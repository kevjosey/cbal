###################################
## TITLE: cbalance_test.R        ##
## PURPOSE: Test run cbalance.R  ##
###################################

library(cbal)
library(survey)

tau <- 20
sig2 <- 5
rho <- 0.5
n <- 200
iter <- 100

# simulate array of data
simDat <- replicate(iter, ks_data(n = n, tau = tau, sig2 = sig2, rho = rho, y_scen = "b", z_scen = "a"))

# fit along with ebal
tauHat_ebal <- tauHat_cbps <- vector(mode = "numeric", length = iter)

for (i in 1:iter) {

  dat <- as.data.frame(simDat[,i])
  treat <- dat$z
  y <- dat$y
  cov_dat <- subset(dat, select = c(x1, x2, x3, x4))
  
  fit_ebal <- cbalance(z ~ x1 + x2 + x3 + x4, data = dat, estimand = "ATT", distance = "entropy")
  fit_cbps <- cbalance(z ~ x1 + x2 + x3 + x4, data = dat, estimand = "ATE", distance = "shifted")
  wts_ebal <- fit_ebal$weights
  wts_cbps <- fit_cbps$weights

  design_ebal <- svydesign(ids = ~ 1, weights = ~ wts_ebal, data = data.frame(y, treat, wts_ebal))
  mod_ebal <- svyglm(y ~ treat, design = design_ebal, family = gaussian)
  tauHat_ebal[i] <- ifelse(fit_ebal$converged, coef(mod_ebal)[2], NA)

  design_cbps <- svydesign(ids = ~ 1, weights = ~ wts_cbps, data = data.frame(y, treat, wts_cbps))
  mod_cbps <- svyglm(y ~ treat, design = design_cbps, family = gaussian)
  tauHat_cbps[i] <- ifelse(fit_cbps$converged, coef(mod_cbps)[2], NA)
  
}

mean(tauHat_ebal)
mean(tauHat_cbps)
