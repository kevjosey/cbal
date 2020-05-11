#' Langrangian Functions for Bregman Distances with Linear Equality Constraints
#'
#' Use \code{lagrange_ent} for the unnormalized relative entropy.
#'
#' Use \code{lagrange_bent} for the binary relative entropy.
#'
#' Use \code{lagrange_sent} for the shifted unnormalized relative entropy.
#'
#' @param coefs vector of Lagrangian multipliers (or the dual vector).
#' @param constraint a matrix that determines the basis of a linear subspace which define the equality constraints 
#' of the optimization problem.
#' @param target the target margins of the linear equality constraints. This vector 
#' should have a length equal to the number of columns in \code{constraint}.
#' @param base_weights a vector of base weights with length equal to the 
#' number of rows in \code{constraint}.
#'
#' @name lagrange
NULL

#' @rdname lagrange
#' @export
lagrange_ent <- function(coefs, constraint, target, base_weights) {
  
  temp <- sum(base_weights*exp(-constraint %*% coefs))
  out <- temp + sum(target * coefs)
  return(out)
  
}

#' @rdname lagrange
#' @export
lagrange_bent <- function(coefs, constraint, target, base_weights) {
  
  weights <- c( base_weights / (base_weights + (1 - base_weights)*exp(constraint %*% coefs)) )
  temp <- sum(weights*log(weights/base_weights) + (1 - weights)*log((1 - weights)/(1 - base_weights)))
  out <- -temp - sum(weights * constraint %*% coefs) + sum(target * coefs)
  return(out)
  
}

#' @rdname lagrange
#' @export
lagrange_sent <- function(coefs, constraint, target, base_weights) {
  
  temp <- sum(constraint %*% coefs - (base_weights - 1)*exp(-constraint %*% coefs))
  out <- -temp + sum(target * coefs)
  return(out)
  
}

#' Horvitz-Thompson Estimating Equations
#' 
#' This function is required for the sandwich variance estimator. It provides the
#' estimating equations for the "meat".
#'
#' @param S the binary vector of sample indicators.
#' @param X the balance functions to be contrained.
#' @param Y the observed responses.
#' @param Z the binary treatment assignment.
#' @param p the estimated balancing weights.
#' @param q the estimated sampling weights.
#' @param base_weights a vector of base weights with length equal to the 
#' number of rows in \code{X}.
#' @param theta target sample means of the balance functions.
#' @param tau causal effect estimate.
#'
#' @name est_eq
NULL

#' @rdname esteq
#' @export
esteq_ATE <- function(X, Y, Z, p, base_weights, tau) {
  
  eq1 <- (2*Z - 1)*p*X
  eq2 <- Z*p*(Y - tau) - (1 - Z)*p*Y
  
  eq <- c(eq1, eq2) 
  return(eq)
  
}

#' @rdname esteq
#' @export
esteq_HTE <- function(X, Y, Z, p, base_weights, tau) {
  
  eq1 <- (2*Z - 1)*p*X
  eq2 <- Z*p*X - base_weights*X
  eq3 <- Z*p*(Y - tau) - (1 - Z)*p*Y
  
  eq <- c(eq1, eq2, eq3) 
  return(eq)
  
}

#' @rdname esteq
#' @export
esteq_transport <- function(S, X, Y, Z, p, base_weights, theta, tau) {
  
  eq1 <- S*(Z*p*X - theta)
  eq2 <- S*((1 - Z)*p*X - theta)
  eq3 <- (1 - S)*(base_weights*X - theta)
  eq4 <- S*p*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}

#' @rdname esteq
#' @export
esteq_fusion <- function(X, Y, Z, S, p, q, base_weights, theta, tau) {
  
  eq1 <- p*q*Z*X - theta
  eq2 <- p*q*(1 - Z)*X - theta
  eq3 <- S*(q*X - theta)
  eq4 <- (1 - S)*(base_weights*X - theta)
  eq5 <- p*q*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4, eq5) 
  return(eq)
  
}

#' Bootstrap Indices
#' 
#' Produces a vector of indices which are utilized in the bootstrap
#' variance estimation functioniality of \code{cestimate()}.
#'
#' @param Z the treatment assignment vector.
#'
#' @export
bootit <- function(Z) {
  
  tidx <- which(Z == 1)
  idx1 <- sample(tidx, size = 1, replace = TRUE)
  
  cidx <- which(Z == 0)
  idx2 <- sample(cidx, size = 1, replace = TRUE)
  
  idx_tmp <- c(idx1, idx2)
  n <- length(Z)
  
  idx3 <- sample.int(n, size = n - length(idx_tmp), replace = TRUE)
  idx <- c(idx_tmp, idx3)
  
}
