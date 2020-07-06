test_that("transport balance", {
  
  gen_data <- function(n, sig2 = 5, y_scen = c("a", "b"), z_scen = c("a", "b"), s_scen = c("a", "b")){
    
    # error variance
    R <- matrix(0, nrow = 2, ncol = 2)
    diag(R) <- 1
    V <- diag(sqrt(sig2), nrow = 2, ncol = 2)
    Sig <- V %*% R %*% V
    
    # covariates
    x1 <- stats::rnorm(n, 0, 1)
    x2 <- stats::rnorm(n, 0, 1)
    x3 <- stats::rnorm(n, 0, 1)
    x4 <- stats::rnorm(n, 0, 1)
    
    # transformed predictors
    u1 <- as.numeric(scale(exp((x1 + x4))))
    u2 <- as.numeric(scale((x1 + x2)^3))
    u3 <- as.numeric(scale((x2 + x3)^2))
    u4 <- as.numeric(scale(log(abs(x3*x4))))
    
    # create matrix
    X <- cbind(int = rep(1, n), x1, x2, x3, x4)
    U <- cbind(int = rep(1, n), u1, u2, u3, u4)
    
    # coefficients
    beta <- c(5, -1, 3, -3, 1)
    alpha <- c(10, -3, -1, 1, 3)
    lambda <- c(0, 0, 1, -0.5, 0.5)
    delta <- c(0, -1, -0.5, 0, 0.5)
    gamma <- c(0, -0.5, -1, -0.5, 1)
    
    # Trial Participation
    if (s_scen == "b") {
      f_X <- 1/(1 + exp( -( U %*% gamma) ) )
    } else { # s_scen == "a"
      f_X <- 1/(1 + exp( -( X %*% gamma) ) )
    } 
    
    s <- rbinom(n, 1, f_X)
    
    # propensity score
    if (z_scen == "b") {
      e_X <- s/(1 + exp( -( U %*% delta) ) ) + (1 - s)/(1 + exp( -( U %*% lambda) ) )
    } else { # z_scen == "a"
      e_X <- s/(1 + exp( -( X %*% delta) ) ) + (1 - s)/(1 + exp( -( X %*% lambda) ) )
    }
    
    z <- rbinom(n, 1, e_X)
    
    if (y_scen == "b") {
      W <- U
    } else { # y_scen == "b"
      W <- X
    }
    
    # outcome mean
    mu_0 <- W %*% beta
    mu_1 <- W %*% alpha
    
    # potential outcomes
    eval <- eigen(Sig, symmetric = TRUE)
    y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
    y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
    y_pot <- y_tmp + cbind(mu_0, mu_1) # include causal effect
    
    # observed outcome
    y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
    
    PATE <- mean(y_pot[s == 0,2] - y_pot[s == 0,1])
    
    # create simulation dataset
    sim <- list(y = y, z = z, X = X, s = s, PATE = PATE)
    
    return(sim)
    
  }
  
  set.seed(42)
  
  # simulate array of data
  simDat <- replicate(100, gen_data(n = 1000, sig2 = 10, s_scen = "a", y_scen = "a", z_scen = "b"), simplify = FALSE)
  
  simResult <- lapply(simDat, function(dat) {
    
    fit_t <- transport(S = dat$s, Z1 = dat$z[dat$s == 1], X = dat$X)
    est_t <- estimate(fit_t, Y1 = dat$y[dat$s == 1])
    fit_f <- fusion(S = dat$s, Z = dat$z, X = dat$X)
    est_f <- estimate(fit_f, Y = dat$y)
    matrix(c(est_t$estimate, est_f$estimate), nrow = 1, dimnames = list(character(0), c("est_t", "est_f")))
    
  })
  
  simResult <- do.call(rbind, simResult)
  
  testthat::expect_equal(round(colMeans(simResult), 3), c(est_t = 3.337, est_f = 3.220))
  
})
