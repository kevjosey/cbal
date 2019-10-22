########################################
## TITLE: hte_test.R                  ##
## PURPOSE: Test run HTE capabilities ##
########################################

library(cbal)

hte_data <- function(n, sig2, rho, y_scen = c("a", "b"), z_scen = c("a", "b")){
  
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
  
  tau <- t(apply(X[z == 1,], 2, mean)) %*% gamma
  
  # potential outcomes
  eval <- eigen(Sig, symmetric = TRUE)
  y_init <- matrix(stats::rnorm(n*2, 0, 1), nrow = n, ncol = 2) # iid potential outcomes
  y_tmp <- t(eval$vectors %*% diag(sqrt(eval$values), nrow = 2) %*% t(y_init)) # SVD
  y_pot <- y_tmp + cbind(mu_0, mu_1) # include causal effect
  
  # observed outcome
  y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
  
  # create simulation dataset
  sim <- list(y = y, z = z, x1 = x1, x2 = x2, x3 = x3, x4 = x4, ps = e_X)
  
  return(sim)
  
}

set.seed(07271989)

tau <- 20
sig2 <- 10
rho <- 0
n <- 1000
iter <- 1000

# simulate array of data
simDat <- replicate(iter, hte_data(n = n, sig2 = sig2, rho = rho, y_scen = "a", z_scen = "a"))

# fit using HTE specification
tau_sent <- vector(mode = "numeric", length = iter)
var_sent <- vector(mode = "numeric", length = iter)
cp_sent <- vector(mode = "numeric", length = iter)

for (i in 1:iter) {
  
  dat <- as.data.frame(simDat[,i])
  Y <- dat$y
  
  fit_sent <- cbalance(z ~ x1 + x2 + x3 + x4, data = dat, distance = "shifted")
  est_sent <- cestimate(fit_sent, Y = Y, method = "sandwich")
  tau_sent[i] <- est_sent$tau
  var_sent[i] <- est_sent$variance
  cp_sent[i] <- tau_sent[i] - sqrt(var_sent[i])*1.96 <= tau & tau_sent[i] + sqrt(var_sent[i])*1.96 >= tau
  
}

mean(tau_sent)
mean(cp_sent)
