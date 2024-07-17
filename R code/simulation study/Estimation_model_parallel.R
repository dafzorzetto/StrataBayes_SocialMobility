
library(parallel)

StrataBayes_Gibbs_parallel_X<- function(data, R=5000, burnin=5000, dim_cluster=10){
  
  return(StrataBayes_Gibbs(X=data[["data"]][["X"]],
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
  
  return(StrataBayes_Gibbs(X=data[["data"]][["X"]][,-c(2:3)],
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

# with P(1)-P(0)
model_setting_1_3 <- mclapply(scenario_1_3, StrataBayes_Gibbs_parallel_X,
                            R=5000, burnin=5000,
                            dim_cluster=10,
                            mc.cores = 5)
model_setting_2_3 <- mclapply(scenario_2_3, StrataBayes_Gibbs_parallel_X,
                              R=5000, burnin=5000,
                              dim_cluster=10,
                              mc.cores = 5)

# with more X
model_setting_1_3_X <- mclapply(scenario_1_3_X, StrataBayes_Gibbs_parallel_X,
                              R=5000, burnin=5000,
                              dim_cluster=10,
                              mc.cores = 5)
model_setting_2_3_X <- mclapply(scenario_2_3_X, StrataBayes_Gibbs_parallel_X,
                              R=5000, burnin=5000,
                              dim_cluster=10,
                              mc.cores = 5)

# with f(X) non linear
model_setting_1_3_2_X <- mclapply(scenario_1_3_2_X, StrataBayes_Gibbs_parallel_X,
                                R=5000, burnin=5000,
                                dim_cluster=10,
                                mc.cores = 5)
model_setting_2_3_2_X <- mclapply(scenario_2_3_2_X, StrataBayes_Gibbs_parallel_X,
                                R=5000, burnin=5000,
                                dim_cluster=10,
                                mc.cores = 5)
model_setting_2_3_2_X_w <- mclapply(scenario_2_3_2_X, StrataBayes_Gibbs_parallel_X_Xw,
                                  R=5000, burnin=5000,
                                  dim_cluster=10,
                                  mc.cores = 5)

# with [P(1)-P(0)]^2
model_setting_1_2 <- mclapply(scenario_1_2, StrataBayes_Gibbs_parallel_X,
                              R=5000, burnin=5000,
                              dim_cluster=10,
                              mc.cores = 5)
model_setting_2_2 <- mclapply(scenario_2_2, StrataBayes_Gibbs_parallel_X,
                              R=5000, burnin=5000,
                              dim_cluster=10,
                              mc.cores = 5)

# with more X
model_setting_1_2_X <- mclapply(scenario_1_2_X, StrataBayes_Gibbs_parallel_X,
                                R=5000, burnin=5000,
                                dim_cluster=10,
                                mc.cores = 5)
model_setting_2_2_X <- mclapply(scenario_2_2_X, StrataBayes_Gibbs_parallel_X,
                                R=5000, burnin=5000,
                                dim_cluster=10,
                                mc.cores = 5)

# with f(X) non linear
model_setting_1_2_2_X <- mclapply(scenario_1_2_2_X, StrataBayes_Gibbs_parallel_X,
                                  R=5000, burnin=5000,
                                  dim_cluster=12,
                                  mc.cores = 5)
model_setting_2_2_2_X <- mclapply(scenario_2_2_2_X, StrataBayes_Gibbs_parallel_X,
                                  R=5000, burnin=5000,
                                  dim_cluster=12,
                                  mc.cores = 5)

#############################################################################

# bias and MSE estimation 

bias_funct<-function(data_sim, model_est){
  
  samples= length(model_est)
  
  ATE_Y <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - data_sim[[s]]$simulated_full$Y_0 -
           (apply(model_est[[s]]$Y1_POST,1, median)-
              apply(model_est[[s]]$Y0_POST,1, median))))
  ATE_P <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_1 - data_sim[[s]]$simulated_full$P_0 -
           (apply(model_est[[s]]$P1_POST,1, median)-
              apply(model_est[[s]]$P0_POST,1, median))))
  bias_P0 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_0 - apply(model_est[[s]]$P0_POST,1, median)))
  bias_P1 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_1 - apply(model_est[[s]]$P1_POST,1, median))) 
  bias_Y0 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_0 - apply(model_est[[s]]$Y0_POST,1, median)))
  bias_Y1 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - apply(model_est[[s]]$Y1_POST,1, median)))
  
  return(list(bias_ATE_Y = ATE_Y, bias_ATE_P = ATE_P, 
              bias_P0 = bias_P0, bias_P1=bias_P1, 
              bias_Y0=bias_Y0, bias_Y1=bias_Y1))
}

bias_1 = bias_funct(data_sim = scenario_1, 
                    model_est=model_setting_1)
bias_2_correct = bias_funct(data_sim = scenario_2, 
                            model_est=model_setting_2_correct)
bias_2_all = bias_funct(data_sim = scenario_2, 
                        model_est=model_setting_2_all)

bias_1_3 = bias_funct(data_sim = scenario_1_3, 
                    model_est=model_setting_1_3)
bias_1_2 = bias_funct(data_sim = scenario_1_2, 
                      model_est=model_setting_1_2)
bias_1_3_X = bias_funct(data_sim = scenario_1_3_X, 
                      model_est=model_setting_1_3_X)
bias_1_2_X = bias_funct(data_sim = scenario_1_2_X, 
                        model_est=model_setting_1_2_X)
bias_1_3_2_X = bias_funct(data_sim = scenario_1_3_2_X, 
                        model_est=model_setting_1_3_2_X)
bias_1_2_2_X = bias_funct(data_sim = scenario_1_2_2_X, 
                        model_est=model_setting_1_2_2_X)

bias_2_3 = bias_funct(data_sim = scenario_2_3, 
                      model_est=model_setting_2_3)
bias_2_2 = bias_funct(data_sim = scenario_2_2, 
                      model_est=model_setting_2_2)
bias_2_3_X = bias_funct(data_sim = scenario_2_3_X, 
                        model_est=model_setting_2_3_X)
bias_2_2_X = bias_funct(data_sim = scenario_2_2_X, 
                        model_est=model_setting_2_2_X)
bias_2_3_2_X = bias_funct(data_sim = scenario_2_3_2_X, 
                          model_est=model_setting_2_3_2_X)
bias_2_2_2_X = bias_funct(data_sim = scenario_2_2_2_X, 
                          model_est=model_setting_2_2_2_X)

sapply(bias_1,function(i) summary(i))
sapply(bias_2_correct,function(i) summary(i))
sapply(bias_2_all,function(i) summary(i))


plot_var_comparison<-function(sim, estimation, n_scenario, V=NULL){
  plot(sim$P_0,apply(estimation$P0_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - P(0)"),
       xlab="simulated", ylab="estimated")
  points(sim$P_0[V==1],apply(estimation$P0_POST,1,median)[V==1], col=4)
  points(sim$P_0[V==2],apply(estimation$P0_POST,1,median)[V==2], col=2)
  points(sim$P_0[V==3],apply(estimation$P0_POST,1,median)[V==3], col=3)
  abline(coef = c(0,1), col="red")
  
  plot(sim$P_1,apply(estimation$P1_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - P(1)"),
       xlab="simulated", ylab="estimated")
  points(sim$P_1[V==1],apply(estimation$P1_POST,1,median)[V==1], col=4)
  points(sim$P_1[V==2],apply(estimation$P1_POST,1,median)[V==2], col=2)
  points(sim$P_1[V==3],apply(estimation$P1_POST,1,median)[V==3], col=3)
  abline(coef = c(0,1), col="red")
  
  plot(sim$Y_0,apply(estimation$Y0_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - Y(0)"),
       xlab="simulated", ylab="estimated")
  points(sim$Y_0[V==1],apply(estimation$Y0_POST,1,median)[V==1], col=4)
  points(sim$Y_0[V==2],apply(estimation$Y0_POST,1,median)[V==2], col=2)
  points(sim$Y_0[V==3],apply(estimation$Y0_POST,1,median)[V==3], col=3)
  abline(coef = c(0,1), col="red")
  
  plot(sim$Y_1,apply(estimation$Y1_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - Y(1)"),
       xlab="simulated", ylab="estimated")
  points(sim$Y_1[V==1],apply(estimation$Y1_POST,1,median)[V==1], col=4)
  points(sim$Y_1[V==2],apply(estimation$Y1_POST,1,median)[V==2], col=2)
  points(sim$Y_1[V==3],apply(estimation$Y1_POST,1,median)[V==3], col=3)
  abline(coef = c(0,1), col="red")
}

boxplot_func <-function(data, name_bias, dim_setttings){
  colour_palette=c("#FAD401","#E09600","#E06100","#E02700","#FA31E5","#B600FA","#5D00F5")
  
  boxplot(bias ~ set, data=data, 
          main=name_bias, ylab="bias", xlab="Scenarios",
          col=colour_palette[1:dim_setttings])
  abline(h=0, col="green")
}

boxplot_allbias<-function(biases){
  
  dim_setttings=length(biases)
  
  bias_ATE_Y_all <- cbind(bias=c(sapply(biases, function(j) j$bias_ATE_Y)), 
                          set=rep(1:dim_setttings, each=samples))
  bias_ATE_P_all <- cbind(bias=c(sapply(biases, function(j) j$bias_ATE_P)), 
                          set=rep(1:dim_setttings, each=samples))
  bias_P0_all <- cbind(bias=c(sapply(biases, function(j) j$bias_P0)), 
                       set=rep(1:dim_setttings, each=samples))
  bias_P1_all <- cbind(bias=c(sapply(biases, function(j) j$bias_P1)), 
                       set=rep(1:dim_setttings, each=samples))
  bias_Y0_all <- cbind(bias=c(sapply(biases, function(j) j$bias_Y0)), 
                       set=rep(1:dim_setttings, each=samples))
  bias_Y1_all <- cbind(bias=c(sapply(biases, function(j) j$bias_Y1)), 
                       set=rep(1:dim_setttings, each=samples))
  
  par(mfrow=c(2,3))
  boxplot_func(data=bias_P0_all, 
               name_bias="bias P0", 
               dim_setttings)
  boxplot_func(data=bias_P1_all, 
               name_bias="bias P1", 
               dim_setttings)
  boxplot_func(data=bias_ATE_P_all, 
               name_bias="bias ATE P", 
               dim_setttings)
  boxplot_func(data=bias_Y0_all, 
               name_bias="bias Y0", 
               dim_setttings)
  boxplot_func(data=bias_Y1_all, 
               name_bias="bias Y1", 
               dim_setttings)
  boxplot_func(data=bias_ATE_Y_all, 
               name_bias="bias ATE Y", 
               dim_setttings)
}


par(mfrow=c(2,2))
sample_number=5
plot_var_comparison(sim=scenario_1[[sample_number]]$simulated_full,
                    estimation=model_setting_1[[sample_number]],
                    n_scenario=1)
plot_var_comparison(sim=scenario_2[[sample_number]]$simulated_full,
                    estimation=model_setting_2_correct[[sample_number]],
                    n_scenario=2)
plot_var_comparison(sim=scenario_2[[sample_number]]$simulated_full,
                    estimation=model_setting_2_all[[sample_number]],
                    n_scenario=2)

boxplot_allbias(biases=list(bias_1, bias_1_3, bias_1_2, 
                            bias_1_3_X, bias_1_3_2_X))

boxplot_allbias(biases=list(bias_2_correct, bias_2_all,
                            bias_2_3, bias_2_2, 
                            bias_2_3_X, bias_2_3_2_X))

##############################################################################

# clusters analysis

cluster_estimated<-function(estimation){
  
  V0_vi_clusters <- sapply(1:samples, function(i) 
    max(estimation[[i]]$estimation_partition_V0$cluster_post_vi))
  V0_binder_clusters <- sapply(1:samples, function(i) 
    max(estimation[[i]]$estimation_partition_V0$cluster_post_binder))
  V1_vi_clusters <- sapply(1:samples, function(i) 
    max(estimation[[i]]$estimation_partition_V1$cluster_post_vi))
  V1_binder_clusters <- sapply(1:samples, function(i) 
    max(estimation[[i]]$estimation_partition_V1$cluster_post_binder))
  
  par(mfrow=c(2,2))
  barplot(table(V0_vi_clusters), main="V(0) - VI")
  barplot(table(V0_binder_clusters), main="V(0) - Binder")
  barplot(table(V1_vi_clusters), main="V(1) - VI")
  barplot(table(V1_binder_clusters), main="V(1) - Binder")
  
  return(list(V0_vi_clusters = V0_vi_clusters, V0_binder_clusters = V0_binder_clusters,
              V1_vi_clusters = V1_vi_clusters, V1_binder_clusters = V1_binder_clusters))
}

cluster_estimated(model_setting_1)
cluster_estimated(model_setting_2_all)
cluster_estimated(model_setting_2_correct)

