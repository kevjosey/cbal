
test_that("homogeneous balance", {

  # simulate scenarios
  ks_data <- function(tau, n, sig2, rho, y_scen = c("a", "b"), z_scen = c("a", "b")) {
    
    # covariates
    x1 <- stats::rnorm(n, 0, 1)
    x2 <- stats::rnorm(n, 0, 1)
    x3 <- stats::rnorm(n, 0, 1)
    x4 <- stats::rnorm(n, 0, 1)
    
    # transformed predictors
    u1 <- as.numeric(scale(exp(x1/2)))
    u2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
    u3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
    u4 <- as.numeric(scale((x2 + x4 + 20)^2))
    
    # treatment probabilities
    if (z_scen == "b")
      e_X <- 1/(1 + exp( -(-u1 + 0.5*u2 - 0.25*u3 - 0.1*u4 ) ) )
    else
      e_X <- 1/(1 + exp( -(-x1 + 0.5*x2 - 0.25*x3 - 0.1*x4 ) ) )
    
    z <- as.integer(stats::runif(n) < e_X)
    
    # error variance
    R <- matrix(rho, nrow = 2, ncol = 2)
    diag(R) <- 1
    V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
    Sig <- V %*% R %*% V
    
    if (y_scen == "b")
      mu <- 210 + 27.4*u1 + 13.7*u2 + 13.7*u3 + 13.7*u4
    else
      mu <- 210 + 27.4*x1 + 13.7*x2 + 13.7*x3 + 13.7*x4
    
    eval <- eigen(Sig, symmetric = TRUE)
    y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
    y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
    y_pot <- y_tmp + cbind(mu, mu + tau) # include causal effect
    
    # observed outcome
    y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
    
    # create simulation dataset
    list(y = y, z = z, x1 = x1, x2 = x2, x3 = x3, x4 = x4, ps = e_X, tau = tau)
  }
  
  set.seed(42)
  
  simDat <- replicate(100, ks_data(n = 1000, tau = 20,  sig2 = 10, rho = 0, y_scen = "a", z_scen = "a"),  simplify = FALSE)
  
  simResult <- lapply(simDat, function(dat) {
      
    fit_bent <- balance_OWATE(Z = dat$z, Y = dat$y, X = model.matrix( ~ x1 + x2 + x3 + x4, data = dat))
    fit_bent$estimate
    
  })
  
  simResult <- do.call(c, simResult)
  
  testthat::expect_equal(round(mean(simResult), 3), expected = 20.018)

})

test_that("heterogeneous balance", {
  
  hte_data <- function(tau, n, sig2, rho, y_scen = c("a", "b"), z_scen = c("a", "b")){
    
    # error variance
    R <- matrix(rho, nrow = 2, ncol = 2)
    diag(R) <- 1
    V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
    Sig <- V %*% R %*% V
    
    # covariates
    x1 <- stats::rnorm(n, 0, 1)
    x2 <- stats::rnorm(n, 0, 1)
    x3 <- stats::rnorm(n, 0, 1)
    x4 <- stats::rnorm(n, 0, 1)
    
    # transformed predictors
    u1 <- as.numeric(scale(exp(x1/2)))
    u2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
    u3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
    u4 <- as.numeric(scale((x2 + x4 + 20)^2))
    
    # effect coefficients
    beta <- c(210, 27.4, 13.7, 13.7, 13.7)
    gamma <- c(20, -13.7, 0, 0, 13.7)
    
    # propensity score
    if (z_scen == "b") {
      e_X <- 1/(1 + exp( -(-u1 + 0.5*u2 - 0.25*u3 - 0.1*u4 ) ) )
    } else { # z_scen == "a"
      e_X <- 1/(1 + exp( -(-x1 + 0.5*x2 - 0.25*x3 - 0.1*x4 ) ) )
    }
    
    z <- rbinom(n, 1, e_X)
    
    if (y_scen == "b") {
      X <- cbind(rep(1, times = n), u1, u2, u3, u4)
    } else { # y_scen == "b"
      X <- cbind(rep(1, times = n), x1, x2, x3, x4)
    }
    
    # outcome mean
    mu_0 <- X%*%beta
    mu_1 <- X%*%(beta + gamma)
    
    # potential outcomes
    eval <- eigen(Sig, symmetric = TRUE)
    y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
    y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
    y_pot <- y_tmp + cbind(mu_0, mu_1) # include causal effect
    
    # observed outcome
    y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
    
    # create simulation dataset
    list(y = y, z = z, x1 = x1, x2 = x2, x3 = x3, x4 = x4, ps = e_X, tau = tau) 
  }
  
  set.seed(42)
  
  # simulate array of data
  simDat <- replicate(100, hte_data(tau = 20, n = 1000, sig2 = 10, rho = 0, y_scen = "a", z_scen = "a"),  simplify = FALSE)
  
  simResult <- lapply(simDat, function(dat) {
    
    fit_sent <- balance_ATE(Z = dat$z, Y = dat$y, X = model.matrix( ~ x1 + x2 + x3 + x4, data = dat))
    fit_sent$estimate
    
  })
  
  simResult <- do.call(c, simResult)
  
  testthat::expect_equal(round(mean(simResult), 3), expected = 20.068)
  
})
