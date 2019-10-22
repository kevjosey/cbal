#' Covariate Balancing Weights via Generalized Projections of Bregman Distances
#'
#' The \code{cbalance()} function solves a convex program with linear equality constraints determined by the
#' estimand (\code{estimand}), the criterion function (\code{distance}), and the sampling weights (\code{base_weights}).
#' The function \code{cbalance.fit()} provides a more direct means to solving the convex program. However, 
#' the constraint matrix and target margins must be determined by the user.
#' 
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted.
#' @param data a \code{data.frame}, list, or environment containing the variables in the model. 
#' @param distance the Bregman distance to be optimized. Can either be "entropy" for the relative entropy,
#' "binary" for the binary relative entropy, or "shifted" for the shifted relative entropy. The distance also determines
#' the causal effect estimand. "shifted" produces balancing weights for estimating the average treatment effect,
#' "entropy" for the average treatment effect of the treated, and "binary" for a constant conditional average treatment effect.
#' @param base_weights a vector of optional sampling weights with length equal to the
#' number of rows in \code{cmat} or \code{X}.
#' @param coefs_init optional initialization points for the dual variables. Defaults to a vector of zeros.
#' @param optim_ctrl a list of arguments that will be passed to \code{optim}.
#' @param ... additional arguments.
#'
#' @references
#'
#' Censor Y, Zenios SA (1998). Parallel Optimization: Theory, Algorithms, and Applications. 1st ed. New York:
#' Oxford University Press.
#'
#' @export
cbalance <- function(formula, data,
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
  
  if (!(distance %in% c("entropy", "binary", "shifted")))
    stop("distance must be either \"entropy\", \"binary\", or \"shifted\"")
  
  if (is.null(base_weights)) { # initialize base_weights
    
    if (distance == "binary") {
      base_weights <- rep(1/2, n)
    } else if (distance == "shifted")
      base_weights <- rep(2, n)
    else { # distance == "entropy"
      base_weights <- rep(1, n)
    }
    
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
  
  # initialize coefs
  if (is.null(coefs_init) & distance != "shifted") {
    coefs_init <- rep(0, times = m) 
  } else if (is.null(coefs_init) & distance == "shifted") {
    coefs_init <- rep(0, times = 2*m)
  } else if (distance == "shifted" & length(coefs_init) != 2*m) {
    stop("length(coefs_init) != 2*ncol(X) required for shifted entropy")
  } else if (distance != "shifted" & length(coefs_init) != m) {
    stop("length(coefs_init) != ncol(X)")
  }
  
  if (distance == "entropy") {
    
    cmat <- as.matrix( (1 - Z)*X )
    target <- c( t(Z*X) %*% base_weights )
    
  } else if (distance == "binary") {
    
    cmat <- as.matrix( (2*Z - 1)*X )
    target <- rep(0, times = m)
    
  } else { # distance == "shifted"
    
    cmat <- cbind(as.matrix((2*Z - 1)*X), as.matrix( Z*X ))
    target <- c(rep(0, times = m), c(t(X) %*% base_weights))
    
  }
  
  converged <- FALSE # initialize convergence indicator
  
  # try direct optimization
  cbal_out <- try( cfit(cmat = cmat,
                        target = target,
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
    
    if (distance == "shifted") {
      coefs <- rep(0, times = 2*m)
    } else { # other distance
      coefs <- rep(0, times = m)
    }
    
    maxit <- optim_ctrl$maxit
    reltol <- optim_ctrl$reltol
    
    for (k in 1:maxit) {
      
      vals <- as.vector(abs(t(cmat) %*% temp - target))
      idx <- which.max(vals)
      u <- as.matrix(cmat[,idx])
      tm <- target[idx]
      
      # cbalance to find coefs
      cbalit <- cfit(cmat = u, target = tm, base_weights = temp,
                     coefs_init = coefs_init[idx], distance = distance,
                     optim_ctrl = optim_ctrl, ...)
      
      weights <- cbalit$weights
      coefs[idx] <- coefs[idx] + cbalit$coefs
      
      # indicates convergence
      if (sum(abs( t(cmat) %*% weights - target )) < reltol) {
        converged <- TRUE
        break 
      }
      
      temp <- weights
      
    }
    
  }
  
  if (!converged)
    warning("model failed to converge")
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              formula = formula,
              data = data,
              distance = distance,
              base_weights = base_weights,
              coefs_init = coefs_init,
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cbalance"
  return(out)
  
}

#' @param cmat a matrix that forms the basis of a linear subspace which define the equality constraints of the
#' convex program.
#' @param target the target margins of the linear equality constraints. This vector
#' should have a length equal to the number of columns in \code{cmat}.
#' 
#' @rdname cbalance
#' @export
cfit <- function(cmat, target,
                 distance = c("entropy", "binary", "shifted"),
                 base_weights = NULL,
                 coefs_init = NULL,
                 optim_ctrl = list(maxit = 500, reltol = 1e-10),
                 ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  if (!(distance %in% c("entropy", "binary", "shifted")))
    stop("distance must be either \"entropy\", \"binary\", or \"shifted\"")
  
  if (distance == "binary") {
    fn <- match.fun(lagrange_bent)
  } else if (distance == "shifted") {
    fn <- match.fun(lagrange_sent)
  } else { # distance == "entropy"
    fn <- match.fun(lagrange_ent)
  }
  
  if (is.null(base_weights)) { # initialize base_weights
    
    if (distance == "binary") {
      base_weights <- rep(1/2, nrow(cmat))
    } else if (distance == "shifted") {
      base_weights <- rep(2, nrow(cmat))
    } else { # distance == "entropy"
      base_weights <- rep(1, nrow(cmat))
    }
    
  } else if (length(base_weights) != nrow(cmat)) { 
    stop("length(base_weights) != sample size")
  }
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = ncol(cmat)) 
  } else if (length(coefs_init) != ncol(cmat)) {
    stop("length(coefs_init) != ncol(cmat)")
  }
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  opt <- stats::optim(coefs_init, fn, method = "BFGS",
                      cmat = cmat,
                      base_weights = base_weights,
                      target = target,
                      control = optim_ctrl, ...)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  
  if (distance == "binary") {
    weights <- c( base_weights / (base_weights + (1 - base_weights)*exp(cmat %*% coefs)) )
  } else if (distance == "shifted") {
    weights <- c( 1 + (base_weights - 1)*exp(-cmat %*% coefs) )
  } else { # distance == "entropy"
    weights <- c( base_weights*exp(-cmat %*% coefs) )
  }
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              cmat = cmat,
              target = target,
              distance = distance,
              base_weights = base_weights, 
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}
