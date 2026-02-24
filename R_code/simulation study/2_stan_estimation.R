############################################################################
#       --- estimation DOUBLE_probitDDP ---
############################################################################

# libraries
library(rstan)
library(parallel)

# Stan code
setwd("../src")
probitDP_double_PS_ps <- stan_model("double_probit_PS_ps.stan")
probitDP_double_med_ps <- stan_model("double_probit_Med_ps.stan")

# information for model estimation:
L <- 8                      # DP truncation point fro both DDP
iter_Stan <- 1800           # total iterations for Stan models
burnin <- 1000

###########################################################################
#                       useful functions
###########################################################################

#  1- initialitation function 
#  2 - prepar data for Stan and estimation

estimation_Stan_PY <- function(data, model){
  
  #propensity score
  p.score <- glm(data$Tr ~ data$X,
                 family = binomial,
                 data = as.data.frame(cbind(data$Tr, data$X)))
  pihat <- predict(p.score, as.data.frame(data$X))
  
  # data_model
  data_stan <- list(P = data$P_obs,
                    Y = data$Y_obs,
                    X = apply(data$X,2, scale), 
                    Tr = data$Tr,
                    ps_hat = pihat,
                    N = length(data$Tr),
                    Q = dim(data$X)[2],
                    L_p = L,
                    L_y = L)
  
  init_function <- function() {
    
    
    list(
      # Stick-breaking parameters (keep small to avoid extreme weights)
      beta0_p = rnorm(data_stan$L_p - 1, -0.5, 0.5),
      beta1_p = matrix(rnorm((data_stan$L_p - 1) * data_stan$Q, 0, 0.2), 
                       data_stan$L_p - 1, data_stan$Q),
      
      beta0_y = rnorm(data_stan$L_y - 1, -0.5, 0.5),
      beta1_y = matrix(rnorm((data_stan$L_y - 1) * data_stan$Q, 0, 0.2), 
                       data_stan$L_y - 1, data_stan$Q),
      beta_p_y = rnorm(data_stan$L_y - 1, 0, 00.2),
      
      # Component parameters for P
      alpha0 = rnorm(data_stan$L_p, mean(data_stan$P), sd(data_stan$P) * 0.6),
      alpha1 = matrix(rnorm(data_stan$L_p * data_stan$Q, 0, 0.2), 
                      data_stan$L_p, data_stan$Q),
      alpha2 = rnorm(data_stan$L_p, 0, 0.2),
      alpha_ps = rnorm(data_stan$L_p, 0, 0.2),
      sigma_p = rgamma(data_stan$L_p, 2, 1),
      
      # Component parameters for Y
      gamma0 = rnorm(data_stan$L_y, mean(data_stan$Y), sd(data_stan$Y) * 0.8),
      gamma1 = matrix(rnorm(data_stan$L_y * data_stan$Q, 0, 0.3), 
                      data_stan$L_y, data_stan$Q),
      gamma2 = rnorm(data_stan$L_y, 0, 0.3),
      gamma3 = rnorm(data_stan$L_y, 0, 0.3),
      gamma4 = rnorm(data_stan$L_y, 0, 0.3),
      gamma_ps = rnorm(data_stan$L_y, 0, 0.3),
      sigma_y = rgamma(data_stan$L_y, 2, 1),
      
      # Hyperparameters for P
      mu0_a = mean(data_stan$P),
      mu1_a = rnorm(data_stan$Q, 0, 0.3),
      mu2_a = rnorm(1, 0, 0.3),
      mu_ps_a = rnorm(1, 0, 0.3),
      sigma0_a = rgamma(1, 2, 1),
      sigma1_a = rgamma(data_stan$Q, 2, 1),
      sigma2_a = rgamma(1, 2, 1),
      sigma_ps_a = rgamma(1, 2, 1),
      
      # Hyperparameters for Y
      mu0_g = mean(data_stan$Y),
      mu1_g = rnorm(data_stan$Q, 0, 0.3),
      mu2_g = rnorm(1, 0, 0.3),
      mu3_g = rnorm(1, 0, 0.3),
      mu4_g = rnorm(1, 0, 0.3),
      mu_ps_g = rnorm(1, 0, 0.3),
      sigma0_g = rgamma(1, 2, 1),
      sigma1_g = rgamma(data_stan$Q, 2, 1),
      sigma2_g = rgamma(1, 2, 1),
      sigma3_g = rgamma(1, 2, 1),
      sigma4_g = rgamma(1, 2, 1),
      sigma_ps_g = rgamma(1, 2, 1)
    )
  }
  
  ## run Stan model
  fit_model <- sampling(model, 
                        data = data_stan,
                        init = init_function,
                        iter = iter_Stan,
                        warmup = burnin,
                        chains = 1)
  
  ## extract posterior draws of E[P|X, Tr] for each X in test set.
  post_P0 = rstan::extract(fit_model, pars= c('cond_mean_p0') )$cond_mean_p0
  post_P1 = rstan::extract(fit_model, pars= c('cond_mean_p1') )$cond_mean_p1
  post_Y0 = rstan::extract(fit_model, pars= c('cond_mean_y0') )$cond_mean_y0
  post_Y1 = rstan::extract(fit_model, pars= c('cond_mean_y1') )$cond_mean_y1
  
  return(list(post_P0 = post_P0, post_P1 = post_P1,
              post_Y0 = post_Y0, post_Y1 = post_Y1))
}


###########################################################################
#                       parallelized code
###########################################################################

est_doubleBNP_s1 <- mclapply(scenario_1, function(d) 
  estimation_Stan_PY(data=d$data, model=probitDP_double_PS_ps), mc.cores = 6)
est_doubleBNP_s2 <- mclapply(scenario_2, function(d)
  estimation_Stan_PY(data=d$data, model=probitDP_double_PS_ps), mc.cores = 6)
est_doubleBNP_s3 <- mclapply(scenario_3, function(d) 
  estimation_Stan_PY(data=d$data, model=probitDP_double_med_ps), mc.cores = 6)
est_doubleBNP_s4 <- mclapply(scenario_4, function(d) 
  estimation_Stan_PY(data=d$data, model=probitDP_double_med_ps), mc.cores = 6)


save(est_doubleBNP_s1, file = "results_doubleBNP_s1.RData")
save(est_doubleBNP_s2, file = "results_doubleBNP_s2.RData")
save(est_doubleBNP_s3, file = "results_doubleBNP_s3.RData")
save(est_doubleBNP_s4, file = "results_doubleBNP_s4.RData")


# print one sample
hist_comparison_1sample <- function(data, est){
  par(mfrow=c(4,2))
  hist(data$P_0, nclass=70 )
  hist(colMeans(est$post_P0), nclass=70 )
  hist(data$P_1, nclass=70 )
  hist(colMeans(est$post_P1), nclass=70 )
  
  hist(data$Y_0, nclass=70 )
  hist(colMeans(est$post_Y0), nclass=70 )
  hist(data$Y_1, nclass=70 )
  hist(colMeans(est$post_Y1), nclass=70 )
  
}

hist_comparison_1sample(data=scenario_1[[1]]$simulated_full, est=est_doubleBNP_s1[[1]])
hist_comparison_1sample(data=scenario_2[[1]]$simulated_full, est=est_doubleBNP_s2[[1]])
hist_comparison_1sample(data=scenario_3[[1]]$simulated_full, est=est_doubleBNP_s3[[1]])
hist_comparison_1sample(data=scenario_4[[1]]$simulated_full, est=est_doubleBNP_s4[[1]])

###########################################################################
#           estimation bias and MSE
###########################################################################


est_resullts <- function(data, est){
  
  P0_est = apply(est$post_P0, 2, mean)
  P1_est = apply(est$post_P1, 2, mean)
  Y0_est = apply(est$post_Y0, 2, mean)
  Y1_est = apply(est$post_Y1, 2, mean)
  
  diff_P0 = data$simulated_full$P_0 - P0_est
  diff_P1 = data$simulated_full$P_1 - P1_est
  diff_Y0 = data$simulated_full$Y_0 - Y0_est
  diff_Y1 = data$simulated_full$Y_1 - Y1_est
  
  bias_ATE_P = mean(diff_P1 - diff_P0)
  bias_ATE_Y = mean(diff_Y1 - diff_Y0)
  
  mse_ATE_P = mean((diff_P1 - diff_P0)^2)
  mse_ATE_Y = mean((diff_Y1 - diff_Y0)^2)
  
  return(list(bias_ATE_P = bias_ATE_P, bias_ATE_Y = bias_ATE_Y,
              mse_ATE_P = mse_ATE_P, mse_ATE_Y = mse_ATE_Y))
}


results_doubleBNP_s1 <- lapply(1:samples, function(x) est_resullts(data = scenario_1[[x]],
                                                                   est=est_doubleBNP_s1[[x]]) )
results_doubleBNP_s2 <- lapply(1:samples, function(x) est_resullts(data = scenario_2[[x]],
                                                                   est=est_doubleBNP_s2[[x]]) )
results_doubleBNP_s3 <- lapply(1:samples, function(x) est_resullts(data = scenario_3[[x]],
                                                                   est=est_doubleBNP_s3[[x]]) )
results_doubleBNP_s4 <- lapply(1:samples, function(x) est_resullts(data = scenario_4[[x]],
                                                                   est=est_doubleBNP_s4[[x]]) )


table_results <- function(results){
  
  bias_P <- sapply(results, function(x) x$bias_ATE_P)
  bias_Y <- sapply(results, function(x) x$bias_ATE_Y)
  mse_P <- sapply(results, function(x) x$mse_ATE_P)
  mse_Y <- sapply(results, function(x) x$mse_ATE_Y)
  
  table <- matrix(c(median(bias_P),IQR(bias_P),
                    median(mse_P),IQR(mse_P),
                    median(bias_Y),IQR(bias_Y),
                    median(mse_Y),IQR(mse_Y)),
                  ncol=2,nrow=4, byrow=TRUE)
  colnames(table)<- c("median", "IQR")
  rownames(table)<- c("bias_P", "MSE_P","bias_Y","MSE_Y")
  
  return(table)
}

table_doubleBNP_s1 <- table_results(results_doubleBNP_s1)
table_doubleBNP_s2 <- table_results(results_doubleBNP_s2)
table_doubleBNP_s3 <- table_results(results_doubleBNP_s3)
table_doubleBNP_s4 <- table_results(results_doubleBNP_s4)
