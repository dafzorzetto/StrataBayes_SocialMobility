#sample_1 <- scenario_1[[1]][["data"]]
#OUT_1 = PostTreatGibbs(X=sample_1[["X"]], 
#                       X_w=sample_1[["X"]], 
#                       Tr=sample_1[["Tr"]], 
#                       P=sample_1[["P_obs"]],
#                       Y=sample_1[["Y_obs"]],
#                       R=5000, burnin=5000)

model_setting_1 <- lapply(1:10, function(s) PostTreatGibbs(
  X=scenario_1[[s]][["data"]][["X"]],
  X_w=scenario_1[[s]][["data"]][["X"]], 
  Tr=scenario_1[[s]][["data"]][["Tr"]], 
  P=scenario_1[[s]][["data"]][["P_obs"]],
  Y=scenario_1[[s]][["data"]][["Y_obs"]],
  R=5000, burnin=5000)) 
model_setting_2_correct <- lapply(1:10, function(s) PostTreatGibbs(
  X=scenario_2[[s]][["data"]][["X"]][,c(1,4:5)],
  X_w=scenario_2[[s]][["data"]][["X"]][,1:3], 
  Tr=scenario_2[[s]][["data"]][["Tr"]], 
  P=scenario_2[[s]][["data"]][["P_obs"]],
  Y=scenario_2[[s]][["data"]][["Y_obs"]],
  R=5000, burnin=5000)) 
model_setting_2_all <- lapply(1:10, function(s) PostTreatGibbs(
  X=scenario_2[[s]][["data"]][["X"]],
  X_w=scenario_2[[s]][["data"]][["X"]], 
  Tr=scenario_2[[s]][["data"]][["Tr"]], 
  P=scenario_2[[s]][["data"]][["P_obs"]],
  Y=scenario_2[[s]][["data"]][["Y_obs"]],
  R=5000, burnin=5000)) 
                         

#library(parallel)
#model_setting_1 <- mclapply(scenario_1, PostTreatGibbs,
#                            R=5000, burnin=5000,
#                            sigma_beta=NULL, mu_beta=NULL, 
#                            gamma_1=1, gamma_2=1,
#                            lambda_eta=NULL, mu_eta=NULL,
#                            gamma_y_1=1, gamma_y_2=1,
#                            mu_eps = NULL, sigma_eps = NULL,
#                            dim_cluster=10,
#                            mc.cores = 5)


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

bias_2_all = bias_funct(data_sim = scenario_2, 
                        model_est=model_setting_2_all)


bias_ATE_1 <- sapply(1:samples, function(s)
                 mean(scenario_1[[s]]$simulated_full$Y_1 - scenario_1[[s]]$simulated_full$Y_0 -
                        (apply(model_setting_1[[s]]$Y1_POST,1, median)-
                        apply(model_setting_1[[s]]$Y0_POST,1, median))))

MSE_ATE_1 <- sapply(1:samples, function(s)
  mean((scenario_1[[s]]$simulated_full$Y_1 - scenario_1[[s]]$simulated_full$Y_0 -
         (apply(model_setting_1[[s]]$Y1_POST,1, median)-
            apply(model_setting_1[[s]]$Y0_POST,1, median)))^2))

bias_P0_1 <- sapply(1:samples, function(s)
  mean(scenario_1[[s]]$simulated_full$P_0 - apply(model_setting_1[[s]]$P0_POST,1, median))
  )

bias_P1_1 <- sapply(1:samples, function(s)
  mean(scenario_1[[s]]$simulated_full$P_1 - apply(model_setting_1[[s]]$P1_POST,1, median))
) 


bias_Y0_1 <- sapply(1:samples, function(s)
  mean(scenario_1[[s]]$simulated_full$Y_0 - apply(model_setting_1[[s]]$Y0_POST,1, median))
)


bias_Y1_1 <- sapply(1:samples, function(s)
  mean(scenario_1[[s]]$simulated_full$Y_1 - apply(model_setting_1[[s]]$Y1_POST,1, median))
)



# bias and MSE estimation SCALE
bias_ATE_1_scale <- sapply(1:samples, function(s)
  mean(scale(scenario_1[[s]]$simulated_full$Y_1) - scale(scenario_1[[s]]$simulated_full$Y_0) -
       (scale(apply(model_setting_1[[s]]$Y1_POST,1, median))-
           scale(apply(model_setting_1[[s]]$Y0_POST,1, median)))))

MSE_ATE_1_scale <- sapply(1:samples, function(s)
  mean((scale(scenario_1[[s]]$simulated_full$Y_1) - scale(scenario_1[[s]]$simulated_full$Y_0) -
          (scale(apply(model_setting_1[[s]]$Y1_POST,1, median))-
              scale(apply(model_setting_1[[s]]$Y0_POST,1, median))))^2))

bias_P0_1_scale <- sapply(1:samples, function(s)
  mean(scale(scenario_1[[s]]$simulated_full$P_0) - 
         scale(apply(model_setting_1[[s]]$P0_POST,1, median)))
)

bias_P1_1_scale <- sapply(1:samples, function(s)
  mean(scale(scenario_1[[s]]$simulated_full$P_1) - 
         scale(apply(model_setting_1[[s]]$P1_POST,1, median)))
) 


bias_Y0_1_scale <- sapply(1:samples, function(s)
  mean(scale(scenario_1[[s]]$simulated_full$Y_0) - scale(
    apply(model_setting_1[[s]]$Y0_POST,1, median)))
)


bias_Y1_1_scale <- sapply(1:samples, function(s)
  mean(scale(scenario_1[[s]]$simulated_full$Y_1) - 
         scale(apply(model_setting_1[[s]]$Y1_POST,1, median)))
)
