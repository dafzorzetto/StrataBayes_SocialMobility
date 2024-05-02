
library(parallel)

StrataBayes_Gibbs_parallel_X<- function(data, R=5000, burnin=5000, dim_cluster=10){
  
  return(StrataBayes_Gibbs_parallel(X=data[["data"]][["X"]],
                                    X_w=data[["data"]][["X"]], 
                                    Tr=data[["data"]][["Tr"]], 
                                    P=data[["data"]][["P_obs"]],
                                    Y=data[["data"]][["Y_obs"]],
                                    sigma_beta=NULL, mu_beta=NULL, 
                                    gamma_1=1, gamma_2=1,
                                    lambda_eta=NULL, mu_eta=NULL,
                                    gamma_y_1=1, gamma_y_2=1,
                                    mu_eps = NULL, sigma_eps = NULL,
                                    R=R, burnin=burnin,
                                    dim_cluster=dim_cluster))
}

StrataBayes_Gibbs_parallel_X_Xw<- function(data, R=5000, burnin=5000, dim_cluster=10){
  
  return(StrataBayes_Gibbs_parallel(X=data[["data"]][["X"]][,c(1,4:5)],
                                    X_w=data[["data"]][["X"]][,1:3], 
                                    Tr=data[["data"]][["Tr"]], 
                                    P=data[["data"]][["P_obs"]],
                                    Y=data[["data"]][["Y_obs"]],
                                    sigma_beta=NULL, mu_beta=NULL, 
                                    gamma_1=1, gamma_2=1,
                                    lambda_eta=NULL, mu_eta=NULL,
                                    gamma_y_1=1, gamma_y_2=1,
                                    mu_eps = NULL, sigma_eps = NULL,
                                    R=R, burnin=burnin,
                                    dim_cluster=dim_cluster))
}

model_setting_1 <- mclapply(scenario_1, StrataBayes_Gibbs_parallel_X,
                            R=5000, burnin=5000,
                            dim_cluster=10,
                            mc.cores = 5)
model_setting_2_correct <- mclapply(scenario_2, StrataBayes_Gibbs_parallel_X,
                                    R=5000, burnin=5000,
                                    dim_cluster=10,
                                    mc.cores = 5)
model_setting_2_all <- mclapply(scenario_2, StrataBayes_Gibbs_parallel_X_Xw,
                                R=5000, burnin=5000,
                                dim_cluster=10,
                                mc.cores = 5)



