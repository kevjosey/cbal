#' Covariate Balancing Weights via Generalized Projections of Bregman Distances
#'
#' The \code{cbalance()} function solves a convex program with linear equality constraints determined by the data, the
#' estimand (\code{estimand}), and the sampling weights (\code{base_weights}).
#' 
#' @param X the balance functions to be contrained.
#' @param Z the binary treatment assignment. 
#' @param estimand the causal estimand to be estimated. Must be either "ATE" for the average treatment effect,
#' "ATT" for the average treatment effect, or "OWATE" for the optimally weighted average treatment effect.
#' In \code{cbalance()}, the estimand also determines the distance function supplied to \code{cfit()}.
#' @param base_weights a vector of optional base weights with length equal to the
#' number of rows in \code{X}.
#' @param coefs_init the optional initial values for the dual variables. Default is a vector of zeros with length 
#' equal to number of columns in \code{X}.
#' @param optim_ctrl a list of arguments that will be passed to \code{optim()}.
#' @param ... additional arguments.
#'
#' @references
#'
#' Censor Y, Zenios SA (1998). Parallel Optimization: Theory, Algorithms, and Applications. 1st ed. New York:
#' Oxford University Press.
#'
#' @rdname cbalance
#' @export
cbalance <- function(X, Y, Z, estimand = c("ATE", "ATT", "OWATE"),
                     base_weights = NULL, coefs_init = NULL,
                     optim_ctrl = list(maxit = 500, reltol = 1e-10), 
                     ...) {
  
  n <- length(Z)
  m <- ncol(X)
  
  # error checks
  if(nlevels(factor(Z)) != 2L)
    stop(paste("nlevels(Z) != 2\nnlevels =", nlevels(factor(Z))))
  
  if (!(estimand %in% c("ATE", "ATT", "OWATE")))
    stop("estimand must be either \"ATE\", \"ATT\", or \"OWATE\"")
  
  if (is.null(base_weights)) { # initialize base_weights
    
    if (estimand == "OWATE") {
      base_weights <- rep(1/2, n)
    } else if (estimand == "ATE")
      base_weights <- rep(2, n)
    else { # distance == "entropy"
      base_weights <- rep(1, n)
    }
    
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
  
  if (estimand == "ATT") {
    
    constraint <- as.matrix( (1 - Z)*X )
    target <- c( t(Z*X) %*% base_weights )
    distance <- "entropy"
    
  } else if (estimand == "OWATE") {
    
    constraint <- as.matrix( (2*Z - 1)*X )
    target <- rep(0, times = m)
    distance <- "binary"
    
  } else { # distance == "shifted"
    
    constraint <- cbind(as.matrix((2*Z - 1)*X), as.matrix( Z*X ))
    target <- c(rep(0, times = m), c(t(X) %*% base_weights))
    distance <- "shifted"
    
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

  converged <- FALSE # initialize convergence indicator
  
  # try direct optimization
  fit_out <- try( cfit(constraint = constraint,
                       target = target,
                       base_weights = base_weights,
                       coefs_init = coefs_init,
                       distance = distance,
                       optim_ctrl = optim_ctrl, ...),
                   silent = TRUE )
  
  if (!inherits(fit_out, "try-error")) {
    
    weights <- fit_out$weights
    converged <- fit_out$converged
    coefs <- fit_out$coefs
    
  } else {
    stop("optimization failed.")
  }
  
  if (!converged)
    warning("model failed to converge")
  
  out <- list(X = X, Z = Z,
              weights = weights,
              coefs = coefs,
              converged = converged,
              constraint = constraint,
              target = target,
              distance = distance,
              base_weights = base_weights,
              coefs_init = coefs_init,
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cbalance"
  return(out)
  
}
