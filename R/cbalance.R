#' Covariate Balancing Weights via Generalized Projections of Bregman Distances
#'
#' The \code{cbalance()} function solves a convex program with linear equality constraints determined by the
#' estimand (\code{estimand}), the criterion function (\code{distance}), and the sampling weights (\code{base_weights}).
#' The function \code{cbalance.fit()} provides a more direct means to solving the convex program. However, 
#' the constraint matrix and target margins must be determined by the user.
#'
#' @references
#'
#' Censor Y, Zenios SA (1998). Parallel Optimization: Theory, Algorithms, and Applications. 1st ed. New York:
#' Oxford University Press.
#'
#' @param formula an object of class formula: a symbolic description of the model to be fitted.
#' @param data a \code{data.frame}, list, or environment containing the variables in the model.
#' @param estimand the assumed causal effect estimand. Can either be "ATE" for the average treatment effect,
#' "ATT" for the average treatment effect of the treated, or "ATC" for the average treatment effect of the controls.
#' @param distance the Bregman distance to be optimized. Can either be "entropy" for the relative entropy,
#' "binary" for the binary relative entropy, or "shifted" for the shifted relative entropy.
#' @param base_weights a vector of optional sampling weights with length equal to the
#' number of rows in \code{data} or \code{constr_mat}.
#' @param coefs_init initial values for the Lagrangian multipliers.
#' @param control a list of arguments that will be passed to \code{optim}.
#' @param ... additional arguments.
#'
#' @export
cbalance <- function(formula,
                     data,
                     estimand = c("ATE", "ATT", "ATC"),
                     distance = c("entropy", "binary", "shifted"),
                     base_weights = NULL,
                     control = list(maxit = 500, reltol = 1e-10), 
                     ...) {
  
  data <- as.data.frame(data)[stats::complete.cases(data),]
  formula <- stats::as.formula(formula, env = environment(data))
  yname <- as.character(formula[[2]])
  
  y <- as.factor(data[,yname])
  z <- ifelse(y == levels(y)[1], 0, 1)
  X <- stats::model.matrix(formula, data = data)
  
  # error checks
  if(nlevels(y) != 2L)
    stop(paste("nlevels(y) != 2\nnlevels = ", nlevels(y)))
  
  if (!(estimand %in% c("ATE", "ATT", "ATC")))
    stop("estimand must be either \"ATE\", \"ATT\", or \"ATC\"")
  
  if (is.null(base_weights)) { # initialize base_weights
    
    if (distance == "binary")
      base_weights <- rep(1/2, nrow(X))
    else if (distance == "shifted")
      base_weights <- rep(2, nrow(X))
    else # distance == "entropy"
      base_weights <- rep(1, nrow(X))
    
  } 
  
  else if (length(base_weights) != nrow(X))
    stop("length(base_weights) != sample size")
  
  if (estimand == "ATT") {
    
    constr_mat <- as.matrix( (1 - z)*X )
    target_margins <- c( t(z*X) %*% base_weights )
    
  } else if (estimand == "ATC") {
    
    constr_mat <- as.matrix( z*X )
    target_margins <- c( t((1 - z)*X) %*% base_weights )
    
  } else { # estimand == "ATE"
    
    constr_mat <- as.matrix( (2*z - 1)*X )
    target_margins <- rep(0, ncol(constr_mat))
    
  }
  
  cbalance.fit(constr_mat = constr_mat,
               target_margins = target_margins,
               base_weights = base_weights,
               distance = distance,
               control = control,
               ...)
  
}

#' @param constr_mat a matrix that forms the basis of a linear subspace which defines the equality 
#' constraints of the convex program.
#' @param target_margins the target margins of the linear equality constraints. This vector
#' should have a length equal to the number of columns in \code{constr_mat}.
#' 
#' @rdname cbalance
#' @export
cbalance.fit <- function(constr_mat,
                         target_margins,
                         distance = c("entropy", "binary", "shifted"),
                         base_weights = NULL,
                         coefs_init = NULL,
                         control = list(maxit = 500, reltol = 1e-10),
                         ...) {

  if (!(distance %in% c("entropy", "binary", "shifted")))
    stop("distance must be either \"entropy\", \"binary\", or \"shifted\"")

  if (distance == "binary")
    fn <- match.fun(lagrange_bent)
  else if (distance == "shifted")
    fn <- match.fun(lagrange_sent)
  else # distance == "entropy"
    fn <- match.fun(lagrange_ent)

  if (is.null(base_weights)) { # initialize base_weights

    if (distance == "binary")
      base_weights <- rep(1/2, nrow(constr_mat))
    else if (distance == "shifted")
      base_weights <- rep(2, nrow(constr_mat))
    else # distance == "entropy"
      base_weights <- rep(1, nrow(constr_mat))

  } else if (length(base_weights) != nrow(constr_mat))
    stop("length(base_weights) != sample size")

  if (is.null(coefs_init))
    coefs_init <- rep(0, ncol(constr_mat)) # initialize coefs
  else if(length(coefs_init) != ncol(constr_mat))
    stop("coefs_init needs to have same length as number of covariates")
  
  chkDots(...)

  opt <- stats::optim(coefs_init, fn, method = "BFGS",
                      constr_mat = constr_mat,
                      base_weights = base_weights,
                      target_margins = target_margins,
                      control = control, ...)

  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par

  if (distance == "binary")
    weights <- c( base_weights / (base_weights + (1 - base_weights)*exp(constr_mat %*% coefs)) )
  else if (distance == "shifted")
    weights <- c( 1 + (base_weights - 1)*exp(-constr_mat %*% coefs) )
  else # distance == "entropy"
    weights <- c( base_weights*exp(-constr_mat %*% coefs) )

  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              constr_mat = constr_mat,
              target_margins = target_margins,
              distance = distance,
              base_weights = base_weights,
              coefs_init = coefs_init)

  if (!converged)
    warning("model failed to converge")

  class(out) <- "cbalance.fit"
  return(out)

}
