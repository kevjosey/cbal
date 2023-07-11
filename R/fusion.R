#' Balancing Weights for Data-Fusion
#'
#' The \code{fusion_ATE()} function finds balancing weights for combining datasets to estimate the target population average treatment effect.
#' This function requires complete individual-level data from both samples whereas the transport function only requires complete data from the
#' study sample and covariate data from the target sample. 
#' 
#' @param S the binary vector of sample indicators.
#' @param X the balance functions to be contrained.
#' @param Y the observed responses.
#' @param Z the binary treatment assignment.
#' @param base_weights a vector of optional base weights with length equal to the
#' number of rows in \code{X}.
#' @param optim_ctrl a list of arguments that will be passed to \code{optim()}.
#' @param ... additional arguments.
#' 
#' @references
#' 
#' Josey KP, Yang F, Ghosh D, Raghavan S (2020). "A Calibration Approach to Transportability and 
#' Data-Fusion with Observational Data." arXiv:2008.06615 [stat].
#'
#' @name fusion
NULL

#' @rdname fusion
#' @export
fusion_ATE <- function(S, X, Y, Z, base_weights = NULL, optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  # error checks
  if(nlevels(factor(Z)) != 2L)
    stop(paste("nlevels(Z) != 2\nnlevels =", nlevels(factor(Z))))
  
  # error checks
  if(nlevels(factor(S)) != 2L)
    stop(paste("nlevels(S) != 2\nnlevels =", nlevels(factor(S))))
  
  n1 <- sum(S)
  n0 <- sum(1 - S)
  n <- n1 + n0
  m <- ncol(X)
  
  X0 <- X[S == 0,]
  X1 <- X[S == 1,]
  Z0 <- Z[S == 0]
  Z1 <- Z[S == 1]
  
  distance <- "entropy"
  coefs_init <- rep(0, times = 4*m) 
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(1, n)
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
  
  theta <- colSums(base_weights[S == 0]*X0)/sum(base_weights[S == 0])
  A <- cbind(as.matrix(S*Z*X), as.matrix(S*(1 - Z)*X), as.matrix((1 - S)*Z*X), as.matrix((1 - S)*(1 - Z)*X))
  b <- c(n1*theta, n1*theta, n0*theta, n0*theta)
  
  # try direct optimization
  fit_out <- try( calibrate(constraint = A,
                       target = b,
                       base_weights = base_weights,
                       coefs_init = coefs_init,
                       distance = distance,
                       optim_ctrl = optim_ctrl),
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
  
  tau <- sum(weights*(2*Z - 1)*Y)/sum(weights*Z)
  
  U <- matrix(0, ncol = 5*m, nrow = 5*m)
  v <- rep(0, times = 5*m + 1)
  meat <- matrix(0, ncol = 5*m + 1, nrow = 5*m + 1)
  
  for (i in 1:n) {
    
    U[1:(4*m),1:(4*m)] <- U[1:(4*m),1:(4*m)] - weights[i] * A[i,] %*% t(A[i,])
    U[1:m,(4*m + 1):(5*m)] <- U[1:m,(4*m + 1):(5*m)] - diag(S[i], m, m)
    U[(m + 1):(2*m),(4*m + 1):(5*m)] <- U[(m + 1):(2*m),(4*m + 1):(5*m)] - diag(S[i], m, m)
    U[(2*m + 1):(3*m),(4*m + 1):(5*m)] <- U[(2*m + 1):(3*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    U[(3*m + 1):(4*m),(4*m + 1):(5*m)] <- U[(3*m + 1):(4*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    U[(4*m + 1):(5*m),(4*m + 1):(5*m)] <- U[(4*m + 1):(5*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    
    v[1:(4*m)] <- v[1:(4*m)] - weights[i]*(2*Z[i] - 1)*(Y[i] - Z[i]*tau)*A[i,]
    v[5*m + 1] <- v[5*m + 1] - weights[i]*Z[i]
    
    meat <- meat +  tcrossprod(esteq_fusion(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                            p = weights[i], q = base_weights[i], 
                                            theta = theta, tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = 5*m + 1, ncol = 5*m + 1)
  invbread[1:(5*m),1:(5*m)] <- U
  invbread[5*m + 1, ] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error"))
    stop("inversion of \"bread\" matrix failed")
  
  sandwich <- bread %*% meat %*% t(bread)
  variance <- sandwich[5*m + 1, 5*m + 1]
  
  out <- list(estimate = tau, variance = variance,
              S = S, X = X, Y = Y, Z = Z, weights = weights,
              coefs = coefs, converged = converged,
              constraint = A, target = b,
              distance = distance, base_weights = base_weights,
              coefs_init = coefs_init, optim_ctrl = optim_ctrl)
  
  class(out) <- "fusion"
  return(out)
  
}
