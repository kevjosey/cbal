#' Generate Kang & Schafer (2007) Scenarios
#'
#' The \code{ks_data()} function generates a dataset according to the experimental scenarios that appear in
#' Kang and Schafer (2007). These scenarios are useful for evauluating different covariate 
#' balancing methods under treatment and/or outcome model mispecification.
#'
#' @param tau the causal effect to be estimated.
#' @param n the sample size.
#' @param sig2 the marginal error variance of the potential outcomes. We assume heteroscedasticity.
#' @param rho the correlation coefficient between the potential outcomes.
#' @param y_scen the outcome scenario. Can be either \code{c("a", "b")}
#' @param z_scen the treatment assignment scenario. Can be either  \code{c("a", "b")}
#'
#' @references
#' 
#' Kang, J. D., & Schafer, J. L. (2007). Demystifying double robustness: A
#' comparison of alternative strategies for estimating a population mean from
#' incomplete data. Statistical Science, 22(4), 523-539.
#'
#' @export
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
    e_X <- 1/(1 + exp( -(-u1 + 0.5*u2 - 0.25*u3 - 0.1*u4) ) )
  else
    e_X <- 1/(1 + exp( -(-x1 + 0.5*x2 - 0.25*x3 - 0.1*x4) ) )

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
  sim <- as.data.frame(cbind(y, z, x1, x2, x3, x4, u1, u2, u3, u4))

  return(sim)

}
