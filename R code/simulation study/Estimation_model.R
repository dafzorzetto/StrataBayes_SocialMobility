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


