#' Covariate Balancing Weights via Generalized Projections of Bregman Distances
#'
#' The \code{cbalance()} function solves a convex program with linear equality constraints determined by the
#' estimand (\code{estimand}), the criterion function (\code{distance}), and the sampling weights (\code{base_weights}).
#' The function \code{cbalance.fit()} provides a more direct means to solving the convex program. However, 
#' the constraint matrix and target margins must be determined by the user.
#' 
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted.
#' @param data a \code{data.frame}, list, or environment containing the variables in the model. 
#' @param estimand the assumed causal effect estimand. Can either be "ATE" for the average treatment effect,
#' "ATT" for the average treatment effect of the treated, or "ATC" for the average treatment effect of the controls.
#' @param distance the Bregman distance to be optimized. Can either be "entropy" for the relative entropy,
#' "binary" for the binary relative entropy, or "shifted" for the shifted relative entropy.
#' @param base_weights a vector of optional sampling weights with length equal to the
#' number of rows in \code{constr_mat} or \code{X}.
#' @param optim_ctrl a list of arguments that will be passed to \code{optim}.
#' @param ... additional arguments.
#'
#' @references
#'
#' Censor Y, Zenios SA (1998). Parallel Optimization: Theory, Algorithms, and Applications. 1st ed. New York:
#' Oxford University Press.
#'
#' @export
cbalance <- function(formula,
                     data,
                     estimand = c("ATE", "ATT", "ATC"),
                     distance = c("entropy", "binary", "shifted"),
                     base_weights = NULL,
                     coefs_init = NULL,
                     optim_ctrl = list(maxit = 500, reltol = 1e-10), 
                     ...) {
  
  data <- as.data.frame(data)[stats::complete.cases(data),]
  formula <- stats::as.formula(formula, env = environment(data))
  nombre <- as.character(formula[[2]])
  
  W <- as.factor(data[,nombre])
  Z <- ifelse(W == levels(W)[1], 0, 1)
  X <- stats::model.matrix(formula, data = data)
  n <- length(Z)
  m <- ncol(X)
  
  # error checks
  if(nlevels(W) != 2L)
    stop(paste("nlevels(Z) != 2\nnlevels =", nlevels(W)))
  
  if (!(estimand %in% c("ATE", "ATT", "ATC")))
    stop("estimand must be either \"ATE\", \"ATT\", or \"ATC\"")
  
  if (!(distance %in% c("entropy", "binary", "shifted")))
    stop("distance must be either \"entropy\", \"binary\", or \"shifted\"")
  
  if (is.null(base_weights)) { # initialize base_weights
    
    if (distance == "binary")
      base_weights <- rep(1/2, n)
    else if (distance == "shifted")
      base_weights <- rep(2, n)
    else # distance == "entropy"
      base_weights <- rep(1, n)
    
  } else if (length(base_weights) != n)
    stop("length(base_weights) != sample size")

  if (is.null(coefs_init))
    coefs_init <- rep(0, times = m) # initialize coefs
  else if (length(coefs_init) != m)
    stop("length(coefs_init) != ncol(X)")
  
  if (estimand == "ATT") {
    
    constr_mat <- as.matrix( (1 - Z)*X )
    target_margins <- c( t(Z*X) %*% base_weights )
    
  } else if (estimand == "ATC") {
    
    constr_mat <- as.matrix( Z*X )
    target_margins <- c( t((1 - Z)*X) %*% base_weights )
    
  } else { # estimand == "ATE"
    
    constr_mat <- as.matrix( (2*Z - 1)*X )
    target_margins <- rep(0, times = m)
    
  }
  
  converged <- FALSE # initialize convergence indicator
  
  # try direct optimization
  cbal_out <- try( cfit(constr_mat = constr_mat,
                        target_margins = target_margins,
                        base_weights = base_weights,
                        coefs_init = coefs_init,
                        distance = distance,
                        optim_ctrl = optim_ctrl, ...),
                   silent = TRUE )
  
  if (!inherits(cbal_out, "try-error")) {
    
    weights <- cbal_out$weights
    converged <- cbal_out$converged
    coefs <- cbal_out$coefs
    
  }
  
  # iterative optimization
  if (!converged | inherits(cbal_out, "try-error")) {
    
    temp <- base_weights
    coefs <- rep(0, times = m)
    maxit <- optim_ctrl$maxit
    reltol <- optim_ctrl$reltol
    
    for (k in 1:maxit) {
      
      vals <- as.vector(abs(t(constr_mat) %*% temp - target_margins))
      idx <- which.max(vals)
      u <- as.matrix(constr_mat[,idx])
      tm <- target_margins[idx]
      
      # cbalance to find coefs
      cbalit <- cfit(constr_mat = u,
                     target_margins = tm,
                     base_weights = temp,
                     coefs_init = coefs_init[idx],
                     distance = distance,
                     optim_ctrl = optim_ctrl,
                     ...)
      
      weights <- cbalit$weights
      coefs[idx] <- coefs[idx] + cbalit$coefs
      
      # indicates convergence
      if (sum(abs( t(constr_mat) %*% weights - target_margins )) < reltol) {
        converged <- TRUE
        break 
      }
      
      temp <- weights
      
    }
    
  }
  
  if (!converged)
    warning("model failed to converge")
  
  
  est <- sum((2*Z -1)*weights*Y)/sum(weights*Z)

  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              Z = Z, X = X,
              estimand = estimand,
              distance = distance,
              base_weights = base_weights,
              coefs_init = coefs_init,
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cbalance"
  return(out)
  
}

#' @param constr_mat a matrix that forms the basis of a linear subspace which define the equality constraints of the
#' convex program.
#' @param target_margins the target margins of the linear equality constraints. This vector
#' should have a length equal to the number of columns in \code{constr_mat}.
#' 
#' @rdname cbalance
#' @export
cfit <- function(constr_mat, 
                 target_margins,
                 distance = c("entropy", "binary", "shifted"),
                 base_weights = NULL,
                 coefs_init = NULL,
                 optim_ctrl = list(maxit = 500, reltol = 1e-10),
                 ...) {
  
  if (!is.matrix(constr_mat))
    stop("constr_mat must be a matrix")
  
  if (!is.vector(target_margins))
    stop("target_margins must be a vector")
  
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
  
  # initialize coefs
  if (is.null(coefs_init))
    coefs_init <- rep(0, times = ncol(constr_mat)) 
  else if (length(coefs_init) != ncol(constr_mat))
    stop("length(coefs_init) != ncol(constr_mat)")
  
  chkDots(...)
  
  opt <- stats::optim(coefs_init, fn, method = "BFGS",
                      constr_mat = constr_mat,
                      base_weights = base_weights,
                      target_margins = target_margins,
                      control = optim_ctrl, ...)
  
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
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}
