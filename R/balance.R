#' Covariate Balancing Weights via Generalized Projections of Bregman Distances
#'
#' The \code{balance_ATE()}, \code{balance_ATT()}, and \code{balance_OWATE()} functions solve a convex program with linear equality constraints determined by the data, the
#' estimand (ATE, ATT, or OWATE), and the sampling weights (\code{base_weights}).
#' 
#' @param X the balance functions to be contrained.
#' @param Y the observed responses.
#' @param Z the binary treatment assignment. 
#' @param base_weights a vector of optional base weights with length equal to the
#' number of rows in \code{X}.
#' @param coefs_init the optional initial values for the dual variables. Default is a vector of zeros with length 
#' equal to number of columns in \code{X}.
#' @param optim_ctrl a list of arguments that will be passed to \code{optim()}.
#' @param ... additional arguments.
#'
#' @references
#' 
#' Josey KP, Juarez-Colunga E, Yang F, Ghosh D (2019). "A Framework for Covariate Balance using Bregman
#' Distances." arXiv:1903.00390 [stat].
#'
#' @name balance
NULL

#' @rdname balance
#' @export
balance_ATE <- function(X, Y, Z, base_weights = NULL, coefs_init = NULL,
                        optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  n <- length(Z)
  m <- ncol(X)
  
  # error checks
  if(nlevels(factor(Z)) != 2L)
    stop(paste("nlevels(Z) != 2\nnlevels =", nlevels(factor(Z))))
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(2, n)
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
    
  constraint <- cbind(as.matrix((2*Z - 1)*X), as.matrix( Z*X ))
  target <- c(rep(0, times = m), c(t(X) %*% base_weights))
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = 2*m)
  } else if (length(coefs_init) != 2*m) {
    stop("length(coefs_init) != 2*ncol(X) -- required for ATE")
  }

  converged <- FALSE # initialize convergence indicator
  
  # try direct optimization
  fit_out <- try( calibrate(constraint = constraint,
                            target = target,
                            base_weights = base_weights,
                            coefs_init = coefs_init,
                            distance = "shifted",
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
  
  # estimate
  tau <- sum((2*Z -1)*weights*Y)/sum(weights*Z)
  dweights <- as.vector( -(base_weights - 1)*exp(-constraint %*% coefs) )
  
  U <- matrix(0, ncol = 2*m, nrow = 2*m)
  v <- rep(0, times = 2*m + 1)
  meat <- matrix(0, ncol = 2*m + 1, nrow = 2*m + 1)
  
  for (i in 1:n) {
    
    U[1:m,1:(2*m)] <- U[1:m,1:(2*m)] + (2*Z[i] - 1) * dweights[i] * (X[i,] %*% t(constraint[i,]))
    U[(m+1):(2*m),1:(2*m)] <- U[(m+1):(2*m),1:(2*m)] + Z[i] * dweights[i] * (X[i,] %*% t(constraint[i,]))
    v[1:(2*m)] <- v[1:(2*m)] + (2*Z[i] - 1) * dweights[i] * (Y[i] - Z[i]*tau) * constraint[i,]
    v[2*m + 1] <- v[2*m + 1] - weights[i]*Z[i]
    meat <- meat + tcrossprod(esteq_HTE(X = X[i,], Y = Y[i], Z = Z[i], p = weights[i], 
                                        q = base_weights[i], tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = 2*m + 1, ncol = 2*m + 1)
  invbread[1:(2*m),1:(2*m)] <- U
  invbread[2*m + 1,] <- v
    
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error"))
    stop("inversion of \"bread\" matrix failed")
  
  sandwich <- bread %*% meat %*% t(bread)
  variance <- sandwich[2*m + 1, 2*m + 1]
  
  out <- list(estimate = tau, variance = variance,
              X = X, Y = Y, Z = Z, weights = weights,
              coefs = coefs, converged = converged,
              constraint = constraint, target = target,
              base_weights = base_weights,
              coefs_init = coefs_init,  optim_ctrl = optim_ctrl)
  
  class(out) <- "balance"
  return(out)
  
}

#' @rdname balance
#' @export
balance_ATT <- function(X, Y, Z, base_weights = NULL, coefs_init = NULL,
                        optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  n <- length(Z)
  m <- ncol(X)
  
  # error checks
  if(nlevels(factor(Z)) != 2L)
    stop(paste("nlevels(Z) != 2\nnlevels =", nlevels(factor(Z))))
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(1, n)
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
  
  constraint <- as.matrix( (1 - Z)*X )
  target <- c( t(Z*X) %*% base_weights )
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = m) 
  } else if (length(coefs_init) != m) {
    stop("length(coefs_init) != ncol(X)")
  }
  
  converged <- FALSE # initialize convergence indicator
  
  # try direct optimization
  fit_out <- try( calibrate(constraint = constraint,
                            target = target,
                            base_weights = base_weights,
                            coefs_init = coefs_init,
                            distance = "entropy",
                            optim_ctrl = optim_ctrl),
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
  
  # estimate
  tau <- sum((2*Z -1)*weights*Y)/sum(weights*Z)
  dweights <- as.vector( -base_weights*exp(-constraint %*% coefs) )
    
    U <- matrix(0, ncol = m, nrow = m)
    v <- rep(0, times = m + 1)
    meat <- matrix(0, ncol = m + 1, nrow = m + 1)
    
  for (i in 1:n) {
    
    U[1:m,1:m] <- U[1:m,1:m] + (2*Z[i] - 1) * dweights[i] * (X[i,] %*% t(constraint[i,]))
    v[1:m] <- v[1:m] + (2*Z[i] - 1) * dweights[i] * (Y[i] - Z[i]*tau) * constraint[i,]
    v[m + 1] <- v[m + 1] - weights[i]*Z[i]
    meat <- meat + tcrossprod(esteq_ATE(X = X[i,], Y = Y[i], Z = Z[i], p = weights[i], tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = m + 1, ncol = m + 1)
  invbread[1:m,1:m] <- U
  invbread[m + 1, ] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error"))
    stop("inversion of \"bread\" matrix failed")
  
  sandwich <- bread %*% meat %*% t(bread)
  variance <- sandwich[m + 1, m + 1]
  
  out <- list(estimate = tau, variance = variance,
              X = X, Y = Y, Z = Z, weights = weights,
              coefs = coefs, converged = converged,
              constraint = constraint, target = target,
              base_weights = base_weights,
              coefs_init = coefs_init,  optim_ctrl = optim_ctrl)
  
  class(out) <- "balance"
  return(out)
  
}

#' @rdname balance
#' @export
balance_OWATE <- function(X, Y, Z, base_weights = NULL, coefs_init = NULL,
                          optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  n <- length(Z)
  m <- ncol(X)
  
  # error checks
  if(nlevels(factor(Z)) != 2L)
    stop(paste("nlevels(Z) != 2\nnlevels =", nlevels(factor(Z))))
  
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(1/2, n)
  } else if (length(base_weights) != n) {
    stop("length(base_weights) != sample size")
  }
    
  constraint <- as.matrix( (2*Z - 1)*X )
  target <- rep(0, times = m)
  
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = m) 
  } else if (length(coefs_init) != m) {
    stop("length(coefs_init) != ncol(X)")
  }
  
  converged <- FALSE # initialize convergence indicator
  
  # try direct optimization
  fit_out <- try( calibrate(constraint = constraint,
                            target = target,
                            base_weights = base_weights,
                            coefs_init = coefs_init,
                            distance = "binary",
                            optim_ctrl = optim_ctrl),
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
  
  # estimate
  tau <- sum((2*Z -1)*weights*Y)/sum(weights*Z)
  dweights <- as.vector( -base_weights*(1 - base_weights)*exp(constraint %*% coefs) / 
                           (base_weights + (1 - base_weights)*exp(constraint %*% coefs))^2 )
    
  U <- matrix(0, ncol = m, nrow = m)
  v <- rep(0, times = m + 1)
  meat <- matrix(0, ncol = m + 1, nrow = m + 1)
  
  for (i in 1:n) {
    
    U[1:m,1:m] <- U[1:m,1:m] + (2*Z[i] - 1) * dweights[i] * (X[i,] %*% t(constraint[i,]))
    v[1:m] <- v[1:m] + (2*Z[i] - 1) * dweights[i] * (Y[i] - Z[i]*tau) * constraint[i,]
    v[m + 1] <- v[m + 1] - weights[i]*Z[i]
    meat <- meat + tcrossprod(esteq_ATE(X = X[i,], Y = Y[i], Z = Z[i], p = weights[i], tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = m + 1, ncol = m + 1)
  invbread[1:m,1:m] <- U
  invbread[m + 1, ] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error"))
    stop("inversion of \"bread\" matrix failed")
  
  sandwich <- bread %*% meat %*% t(bread)
  variance <- sandwich[m + 1, m + 1]
  
  out <- list(estimate = tau, variance = variance,
              X = X, Y = Y, Z = Z, weights = weights,
              coefs = coefs, converged = converged,
              constraint = constraint, target = target,
              base_weights = base_weights,
              coefs_init = coefs_init,  optim_ctrl = optim_ctrl)
  
  class(out) <- "balance"
  return(out)
  
}