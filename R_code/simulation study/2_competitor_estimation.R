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
