#sample_1 <- scenario_1[[1]][["data"]]
#OUT_1 = PostTreatGibbs(X=sample_1[["X"]], 
#                       X_w=sample_1[["X"]], 
#                       Tr=sample_1[["Tr"]], 
#                       P=sample_1[["P_obs"]],
#                       Y=sample_1[["Y_obs"]],
#                       R=5000, burnin=5000)

samples=1

model_setting_1 <- lapply(1:samples, function(s) StrataBayes_Gibbs(
  X=scenario_1[[s]][["data"]][["X"]],
  X_w=scenario_1[[s]][["data"]][["X"]], 
  Tr=scenario_1[[s]][["data"]][["Tr"]], 
  P=scenario_1[[s]][["data"]][["P_obs"]],
  Y=scenario_1[[s]][["data"]][["Y_obs"]],
  R=5000, burnin=5000,
  dim_cluster=10)) 
model_setting_2_correct <- lapply(1:samples, function(s) StrataBayes_Gibbs(
  X=scenario_2[[s]][["data"]][["X"]][,c(1,4:5)],
  X_w=scenario_2[[s]][["data"]][["X"]][,1:3], 
  Tr=scenario_2[[s]][["data"]][["Tr"]], 
  P=scenario_2[[s]][["data"]][["P_obs"]],
  Y=scenario_2[[s]][["data"]][["Y_obs"]],
  R=5000, burnin=5000,
  dim_cluster=10)) 
model_setting_2_all <- lapply(1:samples, function(s) StrataBayes_Gibbs(
  X=scenario_2[[s]][["data"]][["X"]],
  X_w=scenario_2[[s]][["data"]][["X"]], 
  Tr=scenario_2[[s]][["data"]][["Tr"]], 
  P=scenario_2[[s]][["data"]][["P_obs"]],
  Y=scenario_2[[s]][["data"]][["Y_obs"]],
  R=5000, burnin=5000,
  dim_cluster=10)) 
                         




#sapply(1:samples, function(s) table(model_setting_1[[s]]$estimation_partition_V0))
#sapply(1:samples, function(s) table(model_setting_1[[s]]$estimation_partition_V1))


# bias and MSE estimation 

bias_funct<-function(data_sim, model_est){
  
  samples= length(model_est)
  
  ATE <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - data_sim[[s]]$simulated_full$Y_0 -
           (apply(model_est[[s]]$Y1_POST,1, median)-
              apply(model_est[[s]]$Y0_POST,1, median))))
  bias_P0 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_0 - apply(model_est[[s]]$P0_POST,1, median)))
  bias_P1 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_1 - apply(model_est[[s]]$P1_POST,1, median))) 
  bias_Y0 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_0 - apply(model_est[[s]]$Y0_POST,1, median)))
  bias_Y1 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - apply(model_est[[s]]$Y1_POST,1, median)))
    
    return(list(bias_ATE = ATE, bias_P0 = bias_P0, bias_P1=bias_P1, bias_Y0=bias_Y0, bias_Y1=bias_Y1))
}

bias_1 = bias_funct(data_sim = scenario_1, 
                        model_est=model_setting_1)
bias_2_correct = bias_funct(data_sim = scenario_2, 
                        model_est=model_setting_2_correct)
bias_2_all = bias_funct(data_sim = scenario_2, 
                        model_est=model_setting_2_all)

sapply(bias_1,function(i) summary(i))
sapply(bias_2_correct,function(i) summary(i))
sapply(bias_2_all,function(i) summary(i))


plot_var_comparison<-function(sim, estimation, n_scenario){
  plot(sim$P_0,apply(estimation$P0_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - P(0)"),
       xlab="simulated", ylab="estimated")
  abline(coef = c(0,1), col="red")
  
  plot(sim$P_1,apply(estimation$P1_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - P(1)"),
       xlab="simulated", ylab="estimated")
  abline(coef = c(0,1), col="red")
  
  plot(sim$Y_0,apply(estimation$Y0_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - Y(0)"),
       xlab="simulated", ylab="estimated")
  abline(coef = c(0,1), col="red")
  
  plot(sim$P_1,apply(estimation$Y1_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - Y(1)"),
       xlab="simulated", ylab="estimated")
  abline(coef = c(0,1), col="red")
}

hist_var_comparison<-function(sim, estimation, n_scenario){
  
  x_mix=min(sim$P_0,apply(estimation$P0_POST,1,median))
  x_max=max(sim$P_0,apply(estimation$P0_POST,1,median))
  hist(sim$P_0, breaks=100, xlim=c(x_min,x_max),
       main=paste0("scenario ", n_scenario, " - P(0)"),
       xlab="simulated")
  hist(apply(estimation$P0_POST,1,median), 
       breaks=100, xlim=c(x_min,x_max),
       main=paste0("scenario ", n_scenario, " - P(0)"),
       xlab=="estimated")

  
  plot(sim$P_1,apply(estimation$P1_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - P(1)"),
       xlab="simulated", ylab="estimated")
  abline(coef = c(0,1), col="red")
  
  plot(sim$Y_0,apply(estimation$Y0_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - Y(0)"),
       xlab="simulated", ylab="estimated")
  abline(coef = c(0,1), col="red")
  
  plot(sim$P_1,apply(estimation$Y1_POST,1,median), 
       pch=19, cex=0.5, main=paste0("scenario ", n_scenario, " - Y(1)"),
       xlab="simulated", ylab="estimated")
  abline(coef = c(0,1), col="red")
}

par(mfrow=c(2,2))
plot_var_comparison(sim=scenario_1[[1]]$simulated_full,
                    estimation=model_setting_1[[1]],
                    n_scenario=1)
plot_var_comparison(sim=scenario_2[[1]]$simulated_full,
                    estimation=model_setting_2_correct[[1]],
                    n_scenario=2)
plot_var_comparison(sim=scenario_2[[1]]$simulated_full,
                    estimation=model_setting_2_all[[1]],
                    n_scenario=2)

