#' Constrained Convex Optimization of Bregman Distances
#'
#' The \code{cfit()} function solves a convex program with linear equality constraints determined by the 
#' constraint matrix \code{constraint}, the estimand (\code{estimand}), and the sampling weights (\code{base_weights}).
#' The function \code{cfit()} provides a more direct means to solving the convex optimization program. However, 
#' the constraint matrix and target margins must be determined by the user.
#' 
#' @param constraint a matrix that forms the basis of a linear subspace which define the equality constraints of the
#' convex program.
#' @param target the target margins of the linear equality constraints. This vector
#' should have a length equal to the number of columns in \code{constraint}.
#' @param distance the Bregman distance to be optimized. Can either be "entropy" for the relative entropy,
#' "binary" for the binary relative entropy, or "shifted" for the shifted relative entropy. 
#' @param base_weights a vector of optional base weights with length equal to the
#' number of rows in \code{constraint}.
#' @param coefs_init the optional initialization values for the dual variables. Default is a vector of zeros with length 
#' equal to number of columns in \code{X}.
#' @param optim_ctrl a list of arguments that will be passed to \code{optim()}.
#' @param ... additional arguments.
#' 
#' @references
#'
#' Censor Y, Zenios SA (1998). Parallel Optimization: Theory, Algorithms, and Applications. 1st ed. New York:
#' Oxford University Press.
#' 
#' @rdname cfit
#' @export
calibrate <- function(constraint, target, distance = c("entropy", "binary", "shifted"), base_weights = NULL,
                 coefs_init = NULL, optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  if (!is.matrix(constraint))
    stop("constraint must be a matrix")
  
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
      base_weights <- rep(1/2, nrow(constraint))
    } else if (distance == "shifted") {
      base_weights <- rep(2, nrow(constraint))
    } else { # distance == "entropy"
      base_weights <- rep(1, nrow(constraint))
    }
    
  } else if (length(base_weights) != nrow(constraint)) { 
    stop("length(base_weights) != sample size")
  }
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = ncol(constraint)) 
  } else if (length(coefs_init) != ncol(constraint)) {
    stop("length(coefs_init) != ncol(constraint)")
  }
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  opt <- stats::optim(coefs_init, fn, method = "BFGS",
                      constraint = constraint,
                      base_weights = base_weights,
                      target = target,
                      control = optim_ctrl, ...)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  
  if (distance == "binary") {
    weights <- c( base_weights / (base_weights + (1 - base_weights)*exp(constraint %*% coefs)) )
  } else if (distance == "shifted") {
    weights <- c( 1 + (base_weights - 1)*exp(-constraint %*% coefs) )
  } else { # distance == "entropy"
    weights <- c( base_weights*exp(-constraint %*% coefs) )
  }
  
  if (!converged)
    warning("model failed to converge")
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              constraint = constraint,
              target = target,
              distance = distance,
              base_weights = base_weights,
              coefs_init = coefs_init,
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}
