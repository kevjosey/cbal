#' Function for Estimating Causal Effects with \code{cbal}.
#' 
#' This function returns the Horvitz-Thompson estimates and either the sandwich variance
#' estimate or a bootstrapped variance estimate from an object of type \code{cbalance}. 
#' 
#' @param obj object of class "cbalance" or "cboost".
#' @param Y numeric outcome vector. The average effect is usually calculated.
#' @param method method for estimating the variance of the treatment effect estimate.
#' @param boot_iter the number of bootstrap resamples if method = "bootstrap" is selected.
#' @param boot_frac The minimum proportion of entries from each treatment group represented
#' in the bootstrap resample.
#' 
#' @export
cestimate <- function(obj, Y, method = c("sandwich", "bootstrap"), 
                      boot_iter = 1000, boot_frac = 0.1, ...) {
  
  if (!inherits(obj, "cbalance"))
    stop("obj must be of class \"cbalance\"")
  
  if (!(method %in% c("sandwich", "bootstrap")))
    stop("method must be either \"sandwich\" or \"bootstrap\"")
  
  data <- as.data.frame(obj$data)[stats::complete.cases(obj$data),]
  formula <- stats::as.formula(obj$formula, env = environment(data))
  nombre <- as.character(formula[[2]])
  
  W <- as.factor(data[,nombre])
  Z <- ifelse(W == levels(W)[1], 0, 1)
  X <- stats::model.matrix(formula, data = data)
  est <- obj$estimand
  dist <- obj$distance
  
  n <- length(Z)
  m <- ncol(X)
  p <- obj$weights
  q <- obj$base_weights
  coefs <- obj$coefs
  
  tau <- sum((2*Z -1)*p*Y)/sum(p*Z)
  
  if (method == "sandwich") {
    
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
    
    U <- matrix(0, ncol = m, nrow = m)
    v <- rep(0, times = m + 1)
    meat <- matrix(0, ncol = m + 1, nrow = m + 1)
    
    for (i in 1:n) {
      
      U[(1:m),(1:m)] <- U[(1:m),(1:m)] + (2*Z[i] - 1) * dweights[i] * (X[i,] %*% t(A[i,]))
      v[(1:m)] <- v[1:m] + (2*Z[i] - 1) * dweights[i] * (Y[i] - Z[i]*tau) * A[i,]
      v[m + 1] <- v[m + 1] - p[i]*Z[i]
      meat <- meat + tcrossprod(esteq(X = X[i,], Y = Y[i], Z = Z[i], weights = p[i], tau = tau))
      
    }
    
    invbread <- matrix(0, nrow = m + 1, ncol = m + 1)
    invbread[1:m,1:m] <- U
    invbread[m + 1, ] <- v
    
    bread <- solve(invbread)
    sandwich <- bread %*% meat %*% t(bread)
    variance <- sandwich[m + 1, m + 1]
    
  } else if (method == "bootstrap") {
    
    boot_iter <- floor(boot_iter)
    
    if (boot_iter <= 0 | length(boot_iter) != 1)
      stop("boot_iter must be a positive integer")
    
    if (length(boot_frac) != 1 | boot_frac <= 0 | boot_frac > 1)
      stop("boot_frac must be a scalar between 0 and 1")
    
    boot_idx <- replicate(boot_iter, bootit(Z, boot_frac = boot_frac))
    boot_est <- apply( boot_idx, 2, function(idx, obj, Y, Z) {
      
      fit <- cbalance(formula = obj$formula, data = obj$data[idx, ], 
                      estimand = obj$estimand, distance = obj$distance,
                      base_weights = obj$base_weights, coefs_init = obj$coefs_init,
                      obj$optim_ctrl)
      
      r <- fit$weights
      est <- sum((2*Z[idx] - 1)*r*Y[idx])/sum(r*Z[idx])
      return(est)
      
    }, obj = obj, Y = Y, Z = Z )
    
    variance <- var(boot_est)
    
  }
  
  out <- list(tau = tau, variance = variance)
  
}
