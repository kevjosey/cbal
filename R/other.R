#' Langrangian Functions for Bregman Distances with Linear Equality Constraints
#'
#' Use \code{lagrange_ent} for the unnormalized relative entropy.
#'
#' Use \code{lagrange_bent} for the binary relative entropy.
#'
#' Use \code{lagrange_sent} for the shifted unnormalized relative entropy.
#'
#' @param coefs vector of Lagrangian multipliers.
#' @param cmat a matrix that determines the basis of a linear subspace which define the equality constraints 
#' of the optimization problem.
#' @param target the target margins of the linear equality constraints. This vector 
#' should have a length equal to the number of columns in \code{cmat}.
#' @param base_weights a vector of sampling weights with length equal to the 
#' number of rows in \code{cmat}.
#'
#' @name lagrange
NULL

#' @rdname lagrange
#' @export
lagrange_ent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(base_weights*exp(-cmat %*% coefs))
  out <- temp + sum(target * coefs)
  return(out)
  
}

#' @rdname lagrange
#' @export
lagrange_bent <- function(coefs, cmat, target, base_weights) {
  
  weights <- c( base_weights / (base_weights + (1 - base_weights)*exp(cmat %*% coefs)) )
  temp <- sum(weights*log(weights/base_weights) + (1 - weights)*log((1 - weights)/(1 - base_weights)))
  out <- -temp - sum(weights * cmat %*% coefs) + sum(target * coefs)
  return(out)
  
}

#' @rdname lagrange
#' @export
lagrange_sent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(cmat %*% coefs - (base_weights - 1)*exp(-cmat %*% coefs))
  out <- -temp + sum(target * coefs)
  return(out)
  
}

#' Horvitz-Thompson Estimating Equations
#' 
#' This function is required for the sandwich variance estimator. It provides the
#' estimating equations for the "meat".
#'
#' @param X design matrix.
#' @param Y outcome vector.
#' @param Z treatment assignment vector.
#' @param weights estimated balancing weights.
#' @param tau causal effect estimate.
#'
#' @name estimate
NULL

#' @rdname estimate
#' @export
esteq_ATE <- function(X, Y, Z, weights, tau) {
  
  eq1 <- (2*Z - 1)*weights*X
  eq2 <- Z*weights*(Y - tau) - (1 - Z)*weights*Y
  
  eq <- c(eq1, eq2) 
  return(eq)
  
}

#' @rdname estimate
#' @export
esteq_HTE <- function(X, Y, Z, weights, base_weights, tau) {
  
  eq1 <- (2*Z - 1)*weights*X
  eq2 <- Z*weights*X - base_weights*X
  eq3 <- Z*weights*(Y - tau) - (1 - Z)*weights*Y
  
  eq <- c(eq1, eq2, eq3) 
  return(eq)
  
}

#' Bootstrap Indices
#' 
#' Produces a vector of indices which are utilized in the bootstrap
#' variance estimation functioniality of \code{cestimate()}.
#'
#' @param Z treatment assignment vector.
#'
#' @export
bootit <- function(Z) {
  
  tidx <- which(Z == 1)
  idx1 <- sample(tidx, size = 1, replace = TRUE)
  
  cidx <- which(Z == 0)
  idx2 <- sample(cidx, size = 1, replace = TRUE)
  
  idx_tmp <- c(idx1, idx2)
  
  idx3 <- sample.int(n, size = n - length(idx_tmp), replace = TRUE)
  idx <- c(idx_tmp, idx3)
  
}
