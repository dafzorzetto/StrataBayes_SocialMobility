library(doParallel)
library(foreach)

setwd("C:\\Users\\paolo\\OneDrive\\Documents\\Harvard\\code\\results")


# Define the number of simulations and cores for parallel processing
nsim <- 50
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Function to combine results
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# model_setting_2_correct
system.time({
  model_setting_2_correct <- foreach(s = 1:nsim, .multicombine = TRUE) %dopar% {

    StrataBayes_Gibbs(X=scenario_2_3_2_X[[s]][["data"]][["X"]],
                      X_w=scenario_2_3_2_X[[s]][["data"]][["X"]], 
                      Tr=scenario_2_3_2_X[[s]][["data"]][["Tr"]], 
                      P=scenario_2_3_2_X[[s]][["data"]][["P_obs"]],
                      Y=scenario_2_3_2_X[[s]][["data"]][["Y_obs"]],
                      R=5000, burnin=5000,
                      dim_cluster=10  )
  }
})

# Stop the cluster
stopCluster(cl)

save(model_setting_2_correct, file=paste0("model_setting_2_correct_", nsim, "-samples.RData"))

# Define the number of simulations and cores for parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Function to combine results
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# Run the parallel processing
system.time({
  model_setting_2_all <- foreach(s = 1:nsim, .multicombine = TRUE) %dopar% {

    StrataBayes_Gibbs(X=scenario_2_3_2_X[[s]][["data"]][["X"]][,-c(2:3)],
                      X_w=scenario_2_3_2_X[[s]][["data"]][["X"]][,1:3], 
                      Tr=scenario_2_3_2_X[[s]][["data"]][["Tr"]], 
                      P=scenario_2_3_2_X[[s]][["data"]][["P_obs"]],
                      Y=scenario_2_3_2_X[[s]][["data"]][["Y_obs"]],
                      R=5000, burnin=5000,
                      dim_cluster=10  )
  }
})

# Stop the cluster
stopCluster(cl)

save(model_setting_2_all, file=paste0("model_setting_2_all_", nsim, "-samples.RData"))


#model 3 with 14 covariates
# Define the number of simulations and cores for parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Function to combine results
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# Run the parallel processing
system.time({
  model_setting_3_correct <- foreach(s = 1:nsim, .multicombine = TRUE) %dopar% {
    
    StrataBayes_Gibbs(X=scenario_2_3_3_X[[s]][["data"]][["X"]][,-c(2:3)],
                      X_w=scenario_2_3_3_X[[s]][["data"]][["X"]][,1:3], 
                      Tr=scenario_2_3_3_X[[s]][["data"]][["Tr"]], 
                      P=scenario_2_3_3_X[[s]][["data"]][["P_obs"]],
                      Y=scenario_2_3_3_X[[s]][["data"]][["Y_obs"]],
                      R=5000, burnin=5000,
                      dim_cluster=10  )
  }
})

# Stop the cluster
stopCluster(cl)

save(model_setting_3_correct, file=paste0("model_setting_3_correct_", nsim, "-samples.RData"))

#model 3 with 14 covariates
# Define the number of simulations and cores for parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Function to combine results
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# Run the parallel processing
system.time({
  model_setting_3_all <- foreach(s = 1:nsim, .multicombine = TRUE) %dopar% {
    
    StrataBayes_Gibbs(X=scenario_2_3_3_X[[s]][["data"]][["X"]],
                      X_w=scenario_2_3_3_X[[s]][["data"]][["X"]], 
                      Tr=scenario_2_3_3_X[[s]][["data"]][["Tr"]], 
                      P=scenario_2_3_3_X[[s]][["data"]][["P_obs"]],
                      Y=scenario_2_3_3_X[[s]][["data"]][["Y_obs"]],
                      R=5000, burnin=5000,
                      dim_cluster=10  )
  }
})

# Stop the cluster
stopCluster(cl)

save(model_setting_3_all, file=paste0("model_setting_3_all_", nsim, "-samples.RData"))


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
bias_funct_SLM<-function(data_sim, model_est){
  
  samples= length(model_est)
  
  ATE_Y <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - data_sim[[s]]$simulated_full$Y_0 -
           model_est[[s]]$post_Y_1_imp - model_est[[s]]$post_Y_0_imp))
  ATE_P <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_1 - data_sim[[s]]$simulated_full$P_0 -
           model_est[[s]]$post_P_1_imp - model_est[[s]]$post_P_0_imp))
  bias_P0 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_0 - model_est[[s]]$post_P_0_imp))
  bias_P1 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_1 - model_est[[s]]$post_P_1_imp)) 
  bias_Y0 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_0 - model_est[[s]]$post_Y_0_imp))
  bias_Y1 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - model_est[[s]]$post_Y_1_imp))
  
  return(list(bias_ATE_Y = ATE_Y, bias_ATE_P = ATE_P, 
              bias_P0 = bias_P0, bias_P1=bias_P1, 
              bias_Y0=bias_Y0, bias_Y1=bias_Y1))
}


bias_1 = bias_funct(data_sim = scenario_1_3_2_X, 
                    model_est=model_setting_1)
bias_2_correct = bias_funct(data_sim = scenario_2_3_2_X, 
                            model_est=model_setting_2_correct)
bias_2_all = bias_funct(data_sim = scenario_2_3_2_X, 
                        model_est=model_setting_2_all)
bias_3_correct = bias_funct(data_sim = scenario_2_3_3_X, 
                            model_est=model_setting_3_correct)
bias_3_all = bias_funct(data_sim = scenario_2_3_3_X, 
                        model_est=model_setting_3_all)
bias_1_SLM = bias_funct_SLM(data_sim = scenario_1_3_2_X, 
                    model_est=SLM_1)
bias_2_SLM = bias_funct_SLM(data_sim = scenario_2_3_2_X, 
                            model_est=SLM_2)
bias_3_SLM = bias_funct_SLM(data_sim = scenario_2_3_3_X, 
                            model_est=SLM_3)

sapply(bias_1,function(i) summary(i))
sapply(bias_2_correct,function(i) summary(i))
sapply(bias_2_all,function(i) summary(i))
sapply(bias_3_correct,function(i) summary(i))
sapply(bias_3_all,function(i) summary(i))
sapply(bias_1_SLM,function(i) summary(i))
sapply(bias_2_SLM,function(i) summary(i))
sapply(bias_3_SLM,function(i) summary(i))

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

boxplot_allbias(biases=list(bias_1, bias_1_SLM))
boxplot_allbias(biases=list(bias_2_correct, bias_2_all, bias_2_SLM))
boxplot_allbias(biases=list(bias_3_correct, bias_3_all, bias_3_SLM))


