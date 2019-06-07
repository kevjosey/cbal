###################################
## TITLE: cbalance_test.R        ##
## PURPOSE: Test run cbalance.R  ##
###################################

library(cbal)
library(survey)

tau <- 20
sig2 <- 5
rho <- 0
n <- 200
iter <- 1000

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
  sim <- data.frame(y, z, x1, x2, x3, x4, u1, u2, u3, u4)
  
  return(sim)
  
}

set.seed(07271989)

# simulate array of data
simDat <- replicate(iter, ks_data(n = n, tau = tau, sig2 = sig2, rho = rho, y_scen = "a", z_scen = "a"))

# fit along with ebal
tau_ebal <- tau_cbps <- vector(mode = "numeric", length = iter)
var_ebal <- var_cbps <- vector(mode = "numeric", length = iter)
cp_ebal <- cp_cbps <- vector(mode = "numeric", length = iter)

for (i in 1:iter) {
  
  dat <- as.data.frame(simDat[,i])
  Z <- dat$z
  Y <- dat$y
  X <- model.matrix(~ ., data = subset(dat, select = c(x1, x2, x3, x4)))
  
  fit_ebal <- cbalance(z ~ x1 + x2 + x3 + x4, data = dat, estimand = "ATE", distance = "shifted")
  fit_cbps <- cbalance(z ~ x1 + x2 + x3 + x4, data = dat, estimand = "ATE", distance = "binary")
  est_ebal <- cbal_est(fit_ebal, Y = Y) 
  est_cbps <- cbal_est(fit_cbps, Y = Y) 
  
  tau_cbps[i] <- est_cbps$tau
  var_cbps[i] <- est_cbps$var
  cp_cbps[i] <- tau_cbps[i] - sqrt(var_cbps[i])*1.96 <= tau & tau_cbps[i] + sqrt(var_cbps[i])*1.96 >= tau
  
  tau_ebal[i] <- est_ebal$tau
  var_ebal[i] <- est_ebal$var
  cp_ebal[i] <- tau_ebal[i] - sqrt(var_ebal[i])*1.96 <= tau & tau_ebal[i] + sqrt(var_ebal[i])*1.96 >= tau
  
}

mean(tau_ebal)
mean(tau_cbps)
mean(cp_ebal)
mean(cp_cbps)
