#' Langrangian Functions for Bregman Distances with Linear Equality Constraints
#'
#' Use \code{lagrange_ent} for the unnormalized relative entropy.
#'
#' Use \code{lagrange_bent} for the binary relative entropyl
#'
#' Use \code{lagrange_sent} for the shifted unnormalized relative entropyl
#'
#' @param coefs vector of Lagrange multipliers.
#' @param constr_mat a matrix that determines the basis of a linear subspace where the equality constraints of the
#' optimization lie.
#' @param target_margins the target margins of the linear equality constraints. This vector 
#' should have a length equal to the number of columns in \code{constr_mat}.
#' @param base_weights a vector of optional sampling weights with length equal to the 
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
