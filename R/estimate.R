#' Function for Estimating Causal Effects with \code{cbal}.
#' 
#' This function returns the Horvitz-Thompson estimates and the sandwich variance
#' estimate from an object of type \code{cbalance}. 
#' 
#' @param obj object of class "cbalance".
#' @param Y numeric outcome vector.
#' @param method method for estimating the variance of the treatment effect estimate.
#' @param boot_iter the number of bootstrap resamples if method = "bootstrap" is selected.
#' 
#' @export
cestimate <- function(obj, Y, method = c("sandwich", "bootstrap"), boot_iter = 1000, ...) {
  
  if (!inherits(obj, c("cbalance")))
    stop("obj must be of class \"cbalance\"")
  
  if (!(method %in% c("sandwich", "bootstrap")))
    stop("method must be either \"sandwich\" or \"bootstrap\"")
  
  data <- as.data.frame(obj$data)[stats::complete.cases(obj$data),]
  formula <- stats::as.formula(obj$formula, env = environment(data))
  nombre <- as.character(formula[[2]])
  
  W <- as.factor(data[,nombre])
  Z <- ifelse(W == levels(W)[1], 0, 1)
    
  X <- stats::model.matrix(formula, data = data)
  est <- ifelse(obj$distance == "entropy", "ATT", "ATE")
  dist <- obj$distance

  n <- length(Z)
  m <- ncol(X)
  p <- obj$weights
  q <- obj$base_weights
  coefs <- obj$coefs
  
  tau <- sum((2*Z -1)*p*Y)/sum(p*Z)
  
  if (method == "sandwich") {
    
    if (est == "ATT") {
      A <- as.matrix( (1 - Z)*X )
    } else { # est == "ATE"
      A <- as.matrix( (2*Z - 1)*X )
    }

    if (dist == "shifted")
      A <- cbind(as.matrix((2*Z - 1)*X), as.matrix( Z*X ))
      
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
        meat <- meat + tcrossprod(esteq_ATE(X = X[i,], Y = Y[i], Z = Z[i], weights = p[i], tau = tau))
        
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
        meat <- meat + tcrossprod(esteq_HTE(X = X[i,], Y = Y[i], Z = Z[i], weights = p[i], base_weights = q[i], tau = tau))
        
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
  
  out <- list(tau = tau, variance = variance)
  return(out)
  
}
