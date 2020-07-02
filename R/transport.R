#' Balancing Weights for Transportability and Data-Fusion
#'
#' The \code{transport()} and \code{fusion()} functions find balancing weights for transporting estimates. The transport 
#' requires the complete individual-level data (response, treatment assignment, and covariates) for all units in the observed
#' sample but only the individual-level covariate data from the target sample. The fusion function requires complete individual-level
#' data from both samples.
#' 
#' @param S the binary vector of sample indicators.
#' @param X the balance functions to be contrained.
#' @param Z the binary treatment assignment.
#' @param Z1 the binary treatment assignment for S = 1.
#' @param base_weights a vector of optional base weights with length equal to the
#' number of rows in \code{X}.
#' @param coefs_init the optional initialization values for the dual variables. Default is a vector of zeros with length 
#' equal to number of columns in \code{X}.
#' @param optim_ctrl a list of arguments that will be passed to \code{optim()}.
#' @param ... additional arguments.
#' 
#' @references
#' 
#' Josey KP, Berkowitz SA, Ghosh D, Raghavan S (2020). "Transporting Experimental Results with Entropy
#' Balancing." arXiv:2002.07899 [stat].
#'
#' @rdname transport
#' @export
transport <- function(S, X, Z1, base_weights = NULL, coefs_init = NULL,
                      optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  # error checks
  if(nlevels(factor(Z1)) != 2L)
    stop(paste("nlevels(Z1) != 2\nnlevels =", nlevels(factor(Z))))
  
  # error checks
  if(nlevels(factor(S)) != 2L)
    stop(paste("nlevels(S) != 2\nnlevels =", nlevels(factor(S))))
  
  n_1 <- sum(S)
  n_0 <- sum(1 - S)
  n <- n_1 + n_0
  m <- ncol(X)
  X0 <- X[S == 0,]
  X1 <- X[S == 1,]
  Z <- rep(1, times = n)
  Z[S == 1] <- Z1
  distance <- "entropy"
  
  if (is.null(base_weights)) { # initialize base_weights
      base_weights <- rep(1, n)
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
  
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = 2*m) 
  }
  
  theta <- colMeans(base_weights[S == 0]*X0)
  constraint <- cbind(as.matrix(S*Z*X), as.matrix(S*(1 - Z)*X))
  target <- c(n_1*theta, n_1*theta)
  
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
    stop("optimization failed")
  }
  
  if (!converged)
    warning("model failed to converge")
  
  out <- list(S = S, X = X, Z1 = Z1,
              weights = weights,
              coefs = coefs,
              converged = converged,
              constraint = constraint,
              target = target,
              distance = distance,
              base_weights = base_weights,
              coefs_init = coefs_init,
              optim_ctrl = optim_ctrl)
  
  class(out) <- "transport"
  return(out)
  
}

#' @rdname transport
#' @export
fusion <- function(S, X, Z, base_weights = NULL, coefs_init = NULL,
                   optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  # error checks
  if(nlevels(factor(Z)) != 2L)
    stop(paste("nlevels(Z) != 2\nnlevels =", nlevels(factor(Z))))
  
  # error checks
  if(nlevels(factor(S)) != 2L)
    stop(paste("nlevels(S) != 2\nnlevels =", nlevels(factor(S))))
  
  n_1 <- sum(S)
  n_0 <- sum(1 - S)
  n <- n_1 + n_0
  m <- ncol(X)
  
  X0 <- X[S == 0,]
  X1 <- X[S == 1,]
  Z0 <- Z[S == 0]
  Z1 <- Z[S == 1]
  distance <- "entropy"
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(1, n)
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
  
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = 3*m) 
  }
  
  theta <- colMeans(base_weights[S == 0]*X0)
  constraint <- cbind(as.matrix(Z*X), as.matrix( (1 - Z)*X ), as.matrix(S*X))
  target <- c(n*theta, n*theta, n_1*theta)
  
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
    stop("optimization failed")
  }
  
  if (!converged)
    warning("model failed to converge")
  
  out <- list(S = S, X = X, Z = Z,
              weights = weights,
              coefs = coefs,
              converged = converged,
              constraint = constraint,
              target = target,
              distance = distance,
              base_weights = base_weights,
              coefs_init = coefs_init,
              optim_ctrl = optim_ctrl)
  
  class(out) <- "fusion"
  return(out)
  
}
