###################################
## TITLE: balance_test.R        ##
## PURPOSE: Test run balance.R  ##
###################################

library(cbal)

## Homogeneous Treatment effect

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
  
  r_exp <- stats::runif(n)
  z <- ifelse(r_exp < e_X, 1, 0)
  
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
  sim <- data.frame(y = y, z = z, x1 = x1, x2 = x2, x3 = x3, x4 = x4, ps = e_X)
  
  return(sim)
  
}

set.seed(07271989)

tau <- 20
sig2 <- 5
rho <- 0
n <- 1000
iter <- 1000

# simulate array of data
simDat <- replicate(iter, ks_data(n = n, tau = tau, sig2 = sig2, rho = rho, y_scen = "a", z_scen = "a"))

# fit along with sent
tau_bent <- vector(mode = "numeric", length = iter)
var_bent <- vector(mode = "numeric", length = iter)
cp_bent <- vector(mode = "numeric", length = iter)

for (i in 1:iter) {
  
  dat <- as.data.frame(simDat[,i])
  Z <- dat$z
  Y <- dat$y
  X <- model.matrix(~ x1 + x2 + x3 + x4, data = dat)
  
  fit_bent <- balance(Z = Z, X = X, estimand = "OWATE")
  est_bent <- cestimate(fit_bent, Y = Y, method = "sandwich") 
  
  tau_bent[i] <- est_bent$estimate
  var_bent[i] <- est_bent$variance
  cp_bent[i] <- tau_bent[i] - sqrt(var_bent[i])*1.96 <= tau & tau_bent[i] + sqrt(var_bent[i])*1.96 >= tau
  
}

mean(tau_bent)
mean(cp_sent)

## Heterogeneous Treatment Effect

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
  Z <- dat$z
  X <- model.matrix(~ x1 + x2 + x3 + x4, data = dat)
  
  fit_sent <- balance(Z = Z, X = X, estimand = "ATE")
  est_sent <- cestimate(fit_sent, Y = Y, method = "sandwich")
  tau_sent[i] <- est_sent$estimate
  var_sent[i] <- est_sent$variance
  cp_sent[i] <- tau_sent[i] - sqrt(var_sent[i])*1.96 <= tau & tau_sent[i] + sqrt(var_sent[i])*1.96 >= tau
  
}

mean(tau_sent)
mean(cp_sent)
