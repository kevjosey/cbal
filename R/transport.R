#' Balancing Weights for Transporting Observational Results
#'
#' The \code{transport_ATE()} function find balancing weights for transporting estimates. Transport weights 
#' requires the complete individual-level data (response, treatment assignment, and covariates) for all units in the observed
#' sample but only the individual-level covariate data from the target sample.
#' 
#' @param S the binary vector of sample indicators.
#' @param X the balance functions to be contrained.
#' @param Y1 the observed responses for S = 1.
#' @param Z1 the binary treatment assignment for S = 1.
#' @param base_weights a vector of optional base weights with length equal to the
#' number of rows in \code{X}.
#' @param optim_ctrl a list of arguments that will be passed to \code{optim()}.
#' @param ... additional arguments.
#' 
#' @references
#' 
#' Josey KP, Berkowitz SA, Ghosh D, Raghavan S (2020a). "Transporting Experimental Results with Entropy
#' Balancing." arXiv:2002.07899 [stat].
#' 
#' Josey KP, Yang F, Ghosh D, Raghavan S (2020b). "A Calibration Approach to Transportability and 
#' Data-Fusion with Observational Data." arXiv:2008.06615 [stat].
#'
#' @name transport
NULL

#' @rdname transport
#' @export
transport_ATE <- function(S, X, Y1, Z1, base_weights = NULL, optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  # error checks
  if(nlevels(factor(Z1)) != 2L)
    stop(paste("nlevels(Z1) != 2\nnlevels =", nlevels(factor(Z))))
  
  # error checks
  if(nlevels(factor(S)) != 2L)
    stop(paste("nlevels(S) != 2\nnlevels =", nlevels(factor(S))))
  
  n1 <- sum(S)
  n0 <- sum(1 - S)
  n <- n1 + n0
  m <- ncol(X)
  
  X0 <- X[S == 0,]
  X1 <- X[S == 1,]
  Z <- rep(1, times = n)
  Z[S == 1] <- Z1
  
  distance <- "entropy"
  coefs_init <- rep(0, times = 2*m) 
  
  if (is.null(base_weights)) { # initialize base_weights
      base_weights <- rep(1, n)
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
  
  theta <- colSums(base_weights[S == 0]*X0)/sum(base_weights[S == 0])
  A <- cbind(as.matrix(S*(2*Z - 1)*X), as.matrix(S*X))
  b <- c(rep(0,m), n1*theta)
  
  # try direct optimization
  fit_out <- try( calibrate(constraint = A,
                        target = b,
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
  
  Y <- rep(1, times = n)
  Y[S == 1] <- Y1
  Z <- rep(1, times = n)
  Z[S == 1] <- Z1
  
  tau <- sum((S*weights*(2*Z - 1)*Y)/sum(S*Z*weights))
  
  U <- matrix(0, ncol = 3*m, nrow = 3*m)
  v <- rep(0, times = 3*m + 1)
  meat <- matrix(0, ncol = 3*m + 1, nrow = 3*m + 1)
  
  for (i in 1:n) {
    
    U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] - weights[i] * A[i,] %*% t(A[i,])
    U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] - diag(S[i], m, m)
    U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] - diag(1 - S[i], m, m)
    
    v[1:(2*m)] <- v[1:(2*m)] - weights[i]*S[i]*(2*Z[i] - 1)*(Y[i] - Z[i]*tau)*A[i,]
    v[3*m + 1] <- v[3*m + 1] - weights[i]*S[i]*Z[i]
    
    meat <- meat + tcrossprod(esteq_transport(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                               p = weights[i], q = base_weights[i], 
                                               theta = theta, tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = 3*m + 1, ncol = 3*m + 1)
  invbread[1:(3*m),1:(3*m)] <- U
  invbread[3*m + 1, ] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error"))
    stop("inversion of \"bread\" matrix failed")
  
  sandwich <- bread %*% meat %*% t(bread)
  variance <- sandwich[3*m + 1, 3*m + 1]
  
  out <- list(estimate = tau, variance = variance,
              S = S, X = X, Y1 = Y1, Z1 = Z1, weights = weights,
              coefs = coefs, converged = converged,
              constraint = A, target = b,
              distance = distance, base_weights = base_weights,
              coefs_init = coefs_init, optim_ctrl = optim_ctrl)
  
  class(out) <- "transport"
  return(out)
  
}
