#' Function for Estimating Causal Effects with \code{cbal}.
#' 
#' This function returns the Horvitz-Thompson estimates and the sandwich variance
#' estimate from an object of type \code{cbalance}. 
#' 
#' @param obj object of class "cbalance".
#' @param Y numeric outcome vector.
#' 
#' @export
cbal_est <- function(obj, Y, ...) {
  
  if (!inherits(obj, "cbalance"))
    stop("obj must be of class \"cbalance\"")
  
  X <- obj$X
  Z <- obj$Z  
  p <- obj$weights
  q <- obj$base_weights
  coefs <- obj$coefs
  est <- obj$estimand
  dist <- obj$distance
  
  tau <- sum((2*Z -1)*p*Y)/sum(p*Z)
  
  if (est == "ATT")
    A <- (1 - Z)*X
  else if (est == "ATC")
    A <- Z*X
  else # est == "ATE"
    A <- (2*Z - 1)*X
  
  if (dist == "sent")
    dweights <- as.vector( -(q - 1)*exp(-A %*% coefs) )
  else if (dist == "bent")
    dweights <- as.vector( -q*(1 - q)*exp(A %*% coefs) / (q + (1 - q)*exp(A %*% coefs))^2 )
  else # dist == "ent"
    dweights <- as.vector( -q*exp(-A %*% coefs) )
  
  n <- length(Z)
  m <- ncol(A)
  
  U <- matrix(0,ncol = m, nrow = m)
  v <- rep(0, times = m + 1)
  meat <- matrix(0, ncol = m + 1, nrow = m + 1)
  
  for (i in 1:n) {
    
    U[(1:m),(1:m)]<- U[(1:m),(1:m)] + (2*Z[i] - 1) * dweights[i] * (X[i,] %*% t(A[i,]))
    v[(1:m)] <- v[1:m] + (2*Z[i] - 1) * dweights[i] * (Y[i] - Z[i]*tau) * A[i,]
    v[m + 1] <- v[m + 1] - p[i]*Z[i]
    meat <- meat + tcrossprod(est_eq(X = X[i,], Y = Y[i], Z = Z[i], weights = p[i], tau = tau))
    
  }
  
  bread <- matrix(0, nrow = m + 1, ncol = m + 1)
  bread[1:m,1:m] <- U
  bread[m + 1, ] <- v
  
  bread <- bread/n
  meat <- meat/n
  
  bread <- solve(bread)
  var <- (bread %*% meat %*% t(bread))/n
  
  out <- list(tau = tau, var = as.vector(var[m + 1, m + 1]))
  
}
