############################################################################
#       --- estimation DDP model by SLM ---
############################################################################

# upload code 
library(parallel)

# SLM
source("Competitor.R")

#Kim & Zigler
source("BPCF_sample.R")
sourceCpp("code/BPCF/MCMC_main.cpp", rebuild = T, verbose = T)

# information for model estimation:
L <- 8                      # DP truncation point fro both DDP
n_iter <- 2000                # total iterations
burnin <- 1200

###########################################################################
#          useful functions to parallel the estimation
###########################################################################

SLM_Xreg_Gibbs_parallel<- function(data, R=n_iter, burnin=burnin, dim_cluster=L){
  
  return(SLM_Xreg_GIbbs(X=data[["data"]][["X"]],
                        Tr=data[["data"]][["Tr"]], 
                        P=data[["data"]][["P_obs"]],
                        Y=data[["data"]][["Y_obs"]],
                        R=R, burnin=burnin,
                        n_cluster=dim_cluster))
}

# estimation for scenario 1 and scenario 2

SLM_X_1 <- mclapply(scenario_1, SLM_Xreg_Gibbs_parallel,
                    R=n_iter, burnin=burnin,
                    dim_cluster=L,
                    mc.cores = 4)
save(SLM_X_1, file = "SLM_X_1.RData")
SLM_X_2 <- mclapply(scenario_2, SLM_Xreg_Gibbs_parallel,
                    R=n_iter, burnin=burnin,
                    dim_cluster=L,
                    mc.cores = 4)
save(SLM_X_2, file = "SLM_X_2.RData")


# estimation BPCF by Kim and Zigler
BPCF_scenario1 <- lapply(scenario_1, function(x) BPCF_sample(x))
save(BPCF_scenario1, file = "BPCF_1.RData")

BPCF_scenario2 <- lapply(scenario_2, function(x) BPCF_sample(x))
save(BPCF_scenario2, file = "BPCF_2.RData")

BPCF_scenario3 <- lapply(scenario_3, function(x) BPCF_sample(x))
save(BPCF_scenario3, file = "BPCF_3.RData")

BPCF_scenario4 <- lapply(scenario_4, function(x) BPCF_sample(x))
save(BPCF_scenario4, file = "BPCF_4.RData")


###########################################################################
#           estimation bias and MSE
###########################################################################

est_resullts_SLM <- function(data, est){
  
  diff_P0 = data$simulated_full$P_0 - est$post_P_0_imp
  diff_P1 = data$simulated_full$P_1 - est$post_P_1_imp
  diff_Y0 = data$simulated_full$Y_0 - est$post_Y_0_imp
  diff_Y1 = data$simulated_full$Y_1 - est$post_Y_1_imp
  
  bias_ATE_P = mean(diff_P1 - diff_P0)
  bias_ATE_Y = mean(diff_Y1 - diff_Y0)
  
  mse_ATE_P = mean((diff_P1 - diff_P0)^2)
  mse_ATE_Y = mean((diff_Y1 - diff_Y0)^2)
  
  return(list(bias_ATE_P = bias_ATE_P, bias_ATE_Y = bias_ATE_Y,
              mse_ATE_P = mse_ATE_P, mse_ATE_Y = mse_ATE_Y))
}

results_SLM_s1 <- lapply(1:samples, function(x) est_resullts_SLM(data = scenario_1[[x]],
                                                                   est=SLM_X_1[[x]]) )
results_SLM_s2 <- lapply(1:samples, function(x) est_resullts_SLM(data = scenario_2[[x]],
                                                                   est=SLM_X_2[[x]]) )

est_resullts_BPCF <- function(data, est){
  
  diff_P0 = data$simulated_full$P_0 - est$M0_med
  diff_P1 = data$simulated_full$P_1 - est$M1_med
  diff_Y0 = data$simulated_full$Y_0 - est$Y0_med
  diff_Y1 = data$simulated_full$Y_1 - est$Y1_med
  
  bias_ATE_P = mean(diff_P1 - diff_P0)
  bias_ATE_Y = mean(diff_Y1 - diff_Y0)
  
  mse_ATE_P = mean((diff_P1 - diff_P0)^2)
  mse_ATE_Y = mean((diff_Y1 - diff_Y0)^2)
  
  return(list(bias_ATE_P = bias_ATE_P, bias_ATE_Y = bias_ATE_Y,
              mse_ATE_P = mse_ATE_P, mse_ATE_Y = mse_ATE_Y))
}


results_BPCF_s1 <- lapply(1:samples, function(x) est_resullts_BPCF(data = scenario_1[[x]],
                                                                   est=BPCF_scenario1[[x]]) )
results_BPCF_s2 <- lapply(1:samples, function(x) est_resullts_BPCF(data = scenario_2[[x]],
                                                                   est=BPCF_scenario2[[x]]) )
results_BPCF_s3 <- lapply(1:samples, function(x) est_resullts_BPCF(data = scenario_3[[x]],
                                                                   est=BPCF_scenario3[[x]]) )
results_BPCF_s4 <- lapply(1:samples, function(x) est_resullts_BPCF(data = scenario_4[[x]],
                                                                   est=BPCF_scenario4[[x]]) )


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

table_SLM_s1 <- table_results(results_SLM_s1)
table_SLM_s2 <- table_results(results_SLM_s2)

table_BPCF_s1 <- table_results(results_BPCF_s1)
table_BPCF_s2 <- table_results(results_BPCF_s2)
table_BPCF_s3 <- table_results(results_BPCF_s3)
table_BPCF_s4 <- table_results(results_BPCF_s4)
