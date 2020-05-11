#' Function for Estimating Causal Effects using \code{cbal} Objects
#' 
#' This function returns the causal effect estimate and the sandwich variance
#' estimate from an object of type "cbalance", "transport", or "fusion". 
#' 
#' @param obj object of class "cbalance", "transport", or "fusion".
#' @param Y the observed responses.
#' @param Y1 the observed responses for S == 1.
#' @param method method for estimating the variance of the treatment effect estimate.
#' @param boot_iter the number of bootstrap resamples if method = "bootstrap" is selected.
#' @param ... additional arguments.
#' 
#' @rdname estimate
#' @export
cestimate <- function(obj, Y, method = c("sandwich", "bootstrap"), boot_iter = 1000, ...) {
  
  if (!inherits(obj, c("cbalance")))
    stop("obj must be of class \"cbalance\"")
  
  if (!(method %in% c("sandwich", "bootstrap")))
    stop("method must be either \"sandwich\" or \"bootstrap\"")
  
  if (!obj$converged)
    warning("cbalance failed to converge, estimates are not guaranteed to be correct")

  # unpack obj
  Z <- obj$Z
  X <- obj$X
  p <- obj$weights
  q <- obj$base_weights
  coefs <- obj$coefs
  A <- obj$constraint
  dist <- obj$distance
  
  n <- length(Z)
  m <- ncol(X)
  
  tau <- sum((2*Z -1)*p*Y)/sum(p*Z)
  
  if (method == "sandwich") {
      
    if (dist == "shifted") {
      dweights <- as.vector( -(q - 1)*exp(-A %*% coefs) )
    } else if (dist == "binary") {
      dweights <- as.vector( -q*(1 - q)*exp(A %*% coefs) / (q + (1 - q)*exp(A %*% coefs))^2 )
    } else { # dist == "entropy"
      dweights <- as.vector( -q*exp(-A %*% coefs) )
    }
    
    if (dist != "shifted") {
      
      if (dist == "binary") {
        dweights <- as.vector( -q*(1 - q)*exp(A %*% coefs) / (q + (1 - q)*exp(A %*% coefs))^2 )
      } else { # dist == "entropy"
        dweights <- as.vector( -q*exp(-A %*% coefs) )
      }
      
      U <- matrix(0, ncol = m, nrow = m)
      v <- rep(0, times = m + 1)
      meat <- matrix(0, ncol = m + 1, nrow = m + 1)
      
      for (i in 1:n) {
        
        U[1:m,1:m] <- U[1:m,1:m] + (2*Z[i] - 1) * dweights[i] * (X[i,] %*% t(A[i,]))
        v[1:m] <- v[1:m] + (2*Z[i] - 1) * dweights[i] * (Y[i] - Z[i]*tau) * A[i,]
        v[m + 1] <- v[m + 1] - p[i]*Z[i]
        meat <- meat + tcrossprod(esteq_ATE(X = X[i,], Y = Y[i], Z = Z[i], p = p[i], tau = tau))
        
      }
      
      invbread <- matrix(0, nrow = m + 1, ncol = m + 1)
      invbread[1:m,1:m] <- U
      invbread[m + 1, ] <- v
      
      bread <- try(solve(invbread), silent = TRUE)
      
      if (inherits(bread, "try-error"))
        stop("inversion of \"bread\" matrix failed")
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[m + 1, m + 1]
      
    } else {
      
      dweights <- as.vector( -(q - 1)*exp(-A %*% coefs) )
      U <- matrix(0, ncol = 2*m, nrow = 2*m)
      v <- rep(0, times = 2*m + 1)
      meat <- matrix(0, ncol = 2*m + 1, nrow = 2*m + 1)
      
      for (i in 1:n) {
        
        U[1:m,1:(2*m)] <- U[1:m,1:(2*m)] + (2*Z[i] - 1) * dweights[i] * (X[i,] %*% t(A[i,]))
        U[(m+1):(2*m),1:(2*m)] <- U[(m+1):(2*m),1:(2*m)] + Z[i] * dweights[i] * (X[i,] %*% t(A[i,]))
        v[1:(2*m)] <- v[1:(2*m)] + (2*Z[i] - 1) * dweights[i] * (Y[i] - Z[i]*tau) * A[i,]
        v[2*m + 1] <- v[2*m + 1] - p[i]*Z[i]
        meat <- meat + tcrossprod(esteq_HTE(X = X[i,], Y = Y[i], Z = Z[i], p = p[i], base_weights = q[i], tau = tau))
        
      }
      
      invbread <- matrix(0, nrow = 2*m + 1, ncol = 2*m + 1)
      invbread[1:(2*m),1:(2*m)] <- U
      invbread[2*m + 1,] <- v
      
      bread <- try(solve(invbread), silent = TRUE)
      
      if (inherits(bread, "try-error"))
        stop("inversion of \"bread\" matrix failed")
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[2*m + 1, 2*m + 1]
      
    }
    
  } else if (method == "bootstrap") {
    
    boot_iter <- floor(boot_iter)
    
    if (boot_iter <= 0 | length(boot_iter) != 1)
      stop("boot_iter must be a positive integer")
    
    boot_idx <- replicate(boot_iter, bootit(Z))
    boot_est <- apply( boot_idx, 2, function(idx, obj, Y, Z) {
      
      fit <- cbalance(formula = obj$formula, data = obj$data[idx, ], distance = obj$distance,
                      base_weights = obj$base_weights, coefs_init = obj$coefs_init,
                      obj$optim_ctrl)
      
      r <- fit$weights
      est <- sum((2*Z[idx] - 1)*r*Y[idx])/sum(r*Z[idx])
      return(est)
      
    }, obj = obj, Y = Y, Z = Z )
    
    variance <- var(boot_est, na.rm = TRUE)
    
  }
  
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}

#' @rdname estimate
#' @export
testimate <- function(obj, Y1, ...) {
  
  if (!inherits(obj, c("transport")))
    stop("obj must be of class \"transport\"")
  
  S <- obj$S
  X <- obj$X
  Z1 <- obj$Z1
  base_weights <- obj$base_weights
  weights <- obj$weights
  A <- obj$constraint
  coefs <- obj$coefs
  
  n_0 <- sum(1 - S)
  n_1 <- sum(S)
  n <- n_1 + n_0
  m <- ncol(X)
  theta <- obj$target[1:m]/n_1
  
  Y <- rep(1, times = n)
  Y[S == 1] <- Y1
  Z <- rep(1, times = n)
  Z[S == 1] <- Z1
  
  if (is.null(base_weights))
    base_weights <- rep(1, times = length(S))
  
  if (length(base_weights) != length(S))
    stop("base_weights must have the same length as S")
  
  tau <- sum(S*(weights*(2*Z - 1)*Y)/sum(S*Z*weights))
  
  U <- matrix(0, ncol = 3*m, nrow = 3*m)
  v <- rep(0, times = 3*m + 1)
  meat <- matrix(0, ncol = 3*m + 1, nrow = 3*m + 1)
  
  for (i in 1:n) {
    
    U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] - weights[i] * A[i,] %*% t(A[i,])
    
    U[1:m, (2*m + 1):(3*m)] <- U[1:m, (2*m + 1):(3*m)] - diag(S[i], m, m)
    U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] - diag(S[i], m, m)
    U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] - diag((1 - S[i]), m, m)
    
    v[1:(2*m)] <- v[1:(2*m)] - weights[i] * (2*Z[i] - 1) * (Y[i] - Z[i]*tau) * A[i,]
    v[3*m + 1] <- v[3*m + 1] - S[i]*weights[i]*Z[i]
    
    meat <- meat + tcrossprod(esteq_transport(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                              p = weights[i], base_weights = base_weights[i],
                                              theta = theta, tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = 3*m + 1, ncol = 3*m + 1)
  invbread[1:(3*m),1:(3*m)] <- U
  invbread[3*m + 1, ] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error")) {
    
    sandwich <- NA
    variance <- NA
    
  } else {
    
    sandwich <- bread %*% meat %*% t(bread)
    variance <- sandwich[3*m + 1, 3*m + 1]
    
  }
  
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}

#' @rdname estimate
#' @export
festimate <- function(obj, Y, ...) {
  
  if (!inherits(obj, "fusion"))
    stop("obj must be of class \"fusion\"")
  
  S <- obj$S
  X <- obj$X
  Z <- obj$Z
  A <- obj$constraint_p
  R <- obj$constraint_q
  weights <- obj$weights_p
  base_weights <- obj$base_weights
  q <- obj$weights_q
  p <- weights/q
  coefs_1 <- obj$coefs_p
  coefs_2 <- obj$coefs_q
  
  n_0 <- sum(1 - S)
  n_1 <- sum(S)
  n <- n_1 + n_0
  m <- ncol(X)
  
  if (is.null(base_weights))
    base_weights <- rep(1, times = length(S))
  
  if (length(base_weights) != length(S))
    stop("base_weights must have the same length as S")
  
  tau <- sum(weights*(2*Z - 1)*Y)/sum(weights*Z)
  
  dp <- as.vector( -exp(-A %*% coefs_1) )
  dq <- as.vector( -exp(-R %*% coefs_2) )
  U <- matrix(0, ncol = 4*m, nrow = 4*m)
  v <- rep(0, times = 4*m + 1)
  meat <- matrix(0, ncol = 4*m + 1, nrow = 4*m + 1)
  
  for (i in 1:(n_1 + n_0)) {
    
    U[1:(2*m),1:(2*m)] <- U[1:(2*m),1:(2*m)] + q[i]*dp[i] * (A[i,] %*% t(A[i,]))
    
    U[1:(2*m),(2*m + 1):(3*m)] <- U[1:(2*m),(2*m + 1):(3*m)] + 
      dq[i]*p[i]*(A[i,] %*% t(R[i,]))
    U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] + 
      dq[i] * (R[i,] %*% t(R[i,]))
    
    U[1:m,(3*m + 1):(4*m)] <- U[1:m,(3*m + 1):(4*m)] - diag(1, m, m)
    U[(m + 1):(2*m),(3*m + 1):(4*m)] <- U[(m + 1):(2*m),(3*m + 1):(4*m)] - diag(1, m, m)
    U[(2*m + 1):(3*m),(3*m + 1):(4*m)] <- U[(2*m + 1):(3*m),(3*m + 1):(4*m)] - diag(S[i], m, m)
    U[(3*m + 1):(4*m), (3*m + 1):(4*m)] <- U[(3*m + 1):(4*m), (3*m + 1):(4*m)] - diag((1 - S[i]), m, m)
    
    v[1:(2*m)] <- v[1:(2*m)] + q[i]*dp[i] * (2*Z[i] - 1)*(Y[i] - Z[i]*tau)*A[i,]
    v[(2*m+1):(3*m)] <- v[(2*m + 1):(3*m)] + dq[i]*p[i]*(2*Z[i] - 1)*(Y[i] - Z[i]*tau)*R[i,]
    v[4*m + 1] <- v[4*m + 1] - q[i]*p[i]*Z[i]
    
    meat <- meat + tcrossprod(esteq_fusion(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                           p = p[i], q = q[i], base_weights = base_weights[i],
                                           theta = theta, tau = tau))
    
    
  }
  
  invbread <- matrix(0, nrow = 4*m + 1, ncol = 4*m + 1)
  invbread[1:(4*m),1:(4*m)] <- U
  invbread[4*m + 1,] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error"))
    variance <- NA
  
  else {
    
    sandwich <- bread %*% meat %*% t(bread)
    variance <- sandwich[4*m + 1, 4*m + 1]
    
  }
  
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}
