
library(Rcpp)
library(MCMCpack)
library(rootSolve)

sourceCpp("MCMC_main.cpp", rebuild = T, verbose = T)


BPCF_sample <- function(data_sample) {
  set.seed(111)
  n.iter <- 800
  
  # Data 
  Xpred <- data_sample$data$X
  Y_trt <- data_sample$data$Tr
  M_out <- data_sample$data$P_obs
  Y_out <- data_sample$data$Y_obs

  
  PS.fit <- glm(Y_trt~Xpred, family=binomial())
  PS <- predict(PS.fit, type="response")
  
  ### Initial Setup (priors, initial values and hyper-parameters)
  p.grow <- 0.28            # Prob. of GROW
  p.prune <- 0.28           # Prob. of PRUNE
  p.change <- 0.44          # Prob. of CHANGE
  m <- 40                  # Num. of Trees: default setting 100
  
  sigma2_m <- var(M_out)       # Initial value of SD^2
  sigma2_y <- var(Y_out) 
  
  nu <- 3                   # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
  f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_y)
  lambda_y <- rootSolve::uniroot.all(f, c(0.1^5,10))
  f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_m)
  lambda_m <- rootSolve::uniroot.all(f, c(0.1^5,10))
  
  sigma2 <- 1
  f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
  lambda <- rootSolve::uniroot.all(f, c(0.1^5,10))
  
  alpha <- 0.95             # alpha (1+depth)^{-beta} where depth=0,1,2,...
  beta <- 2                 # default setting (alpha, beta) = (0.95, 2)
  
  f <- function(scale) qcauchy(0.75, 0, scale) - 2*sd(M_out-mean(M_out))   # first
  sigma_mu_m_mu_sigma <- uniroot.all(f, c(0.1^5, 100))
  
  f <- function(sd) qnorm(0.75, 0, sd) - sd(M_out-mean(M_out))             # second
  sigma_mu_m_tau_sigma <- uniroot.all(f, c(0.1^5, 100))
  
  f <- function(scale) qcauchy(0.75, 0, scale) - 2*sd(Y_out-mean(Y_out))   # first
  sigma_mu_y_mu_sigma <- uniroot.all(f, c(0.1^5, 100))
  
  f <- function(sd) qnorm(0.75, 0, sd) - sd(Y_out-mean(Y_out))             # second
  sigma_mu_y_tau_sigma <- uniroot.all(f, c(0.1^5, 100))
  
  
  rcpp = MCMC(Xpred, Y_trt, M_out, Y_out, as.numeric(PS),
              p.grow, p.prune, p.change,
              m, 50, m, 50, nu, lambda_m, lambda_y, alpha, beta,
              n.iter,
              sigma_mu_m_tau_sigma, sigma_mu_m_mu_sigma, sigma_mu_y_tau_sigma, sigma_mu_y_mu_sigma)
  
  out <- list(
    S_med  = apply(rcpp$predicted_S,  1, median, na.rm = TRUE),
    Y1_med = apply(rcpp$predicted_Y1, 1, median, na.rm = TRUE),
    Y0_med = apply(rcpp$predicted_Y0, 1, median, na.rm = TRUE),
    M1_med = apply(rcpp$predicted_M1, 1, median, na.rm = TRUE),
    M0_med = apply(rcpp$predicted_M0, 1, median, na.rm = TRUE)
  )
  
  print(paste0("sample ", gg))
  gg = gg + 1
  
  return(out)
}
