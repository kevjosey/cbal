#' Langrangian Functions for Bregman Distances with Linear Equality Constraints
#'
#' Use \code{lagrange_ent} for the unnormalized relative entropy.
#'
#' Use \code{lagrange_bent} for the binary relative entropy.
#'
#' Use \code{lagrange_sent} for the shifted unnormalized relative entropy.
#'
#' @param coefs vector of Lagrangian multipliers.
#' @param constr_mat a matrix that determines the basis of a linear subspace which define the equality constraints 
#' of the optimization problem.
#' @param target_margins the target margins of the linear equality constraints. This vector 
#' should have a length equal to the number of columns in \code{constr_mat}.
#' @param base_weights a vector of sampling weights with length equal to the 
#' number of rows in \code{constr_mat}.
#'
#' @name lagrange
NULL

#' @rdname lagrange
#' @export
lagrange_ent <- function(coefs, constr_mat, target_margins, base_weights) {

  temp <- sum(base_weights*exp(-constr_mat %*% coefs))
  out <- temp + sum(target_margins * coefs)
  return(out)

}

#' @rdname lagrange
#' @export
lagrange_bent <- function(coefs, constr_mat, target_margins, base_weights) {

  weights <- c( base_weights / (base_weights + (1 - base_weights)*exp(constr_mat %*% coefs)) )
  temp <- sum(weights*log(weights/base_weights) + (1 - weights)*log((1 - weights)/(1 - base_weights)))
  out <- -temp - sum(weights * constr_mat %*% coefs) + sum(target_margins * coefs)
  return(out)

}

#' @rdname lagrange
#' @export
lagrange_sent <- function(coefs, constr_mat, target_margins, base_weights) {

  temp <- sum(constr_mat %*% coefs - (base_weights - 1)*exp(-constr_mat %*% coefs))
  out <- -temp + sum(target_margins * coefs)
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
est_eq <- function(X, Y, Z, weights, tau) {
  
  eq1 <- (2*Z - 1)*weights*X
  eq2 <- Z*weights*(Y - tau) - (1 - Z)*weights*Y
  
  eq <- c(eq1, eq2) 
  return(eq)
  
}
