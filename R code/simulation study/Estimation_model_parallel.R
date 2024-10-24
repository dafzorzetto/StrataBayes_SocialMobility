
library(parallel)
source("./src/Gibbs_Sampler_240726.R")
source("./src/competitor_SLM.R")
source("./src/competitor_SLM_Xreg.R")

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

SLM_Gibbs_parallel<- function(data, R=5000, burnin=5000, dim_cluster=10){
  
  return(SLM_GIbbs(X=data[["data"]][["X"]],
                   Tr=data[["data"]][["Tr"]], 
                   P=data[["data"]][["P_obs"]],
                   Y=data[["data"]][["Y_obs"]],
                   R=R+burnin, burnin=burnin,
                   n_cluster=dim_cluster))
}

SLM_Xreg_Gibbs_parallel<- function(data, R=5000, burnin=5000, dim_cluster=10){
  
  return(SLM_Xreg_GIbbs(X=data[["data"]][["X"]],
                        Tr=data[["data"]][["Tr"]], 
                        P=data[["data"]][["P_obs"]],
                        Y=data[["data"]][["Y_obs"]],
                        R=R+burnin, burnin=burnin,
                        n_cluster=dim_cluster))
}

model_setting_1 <- mclapply(scenario_1, StrataBayes_Gibbs_parallel_X,
                            R=2000, burnin=4000,
                            dim_cluster=10,
                            mc.cores = 6)
model_setting_2_a <- mclapply(scenario_2, StrataBayes_Gibbs_parallel_X,
                              R=2000, burnin=4000,
                              dim_cluster=10,
                              mc.cores = 6)
model_setting_2_b <- mclapply(scenario_2, StrataBayes_Gibbs_parallel_X_Xw,
                              R=2000, burnin=4000,
                              dim_cluster=10,
                              mc.cores = 6)
model_setting_3_a <- mclapply(scenario_3, StrataBayes_Gibbs_parallel_X,
                              R=2000, burnin=4000,
                              dim_cluster=10,
                              mc.cores = 6)
model_setting_3_b <- mclapply(scenario_3, StrataBayes_Gibbs_parallel_X_Xw,
                              R=2000, burnin=4000,
                              dim_cluster=10,
                              mc.cores = 6)



prova<-list(scenario_1[[1]],scenario_1[[2]],scenario_1[[3]], scenario_1[[4]], 
            scenario_1[[5]], scenario_1[[6]],scenario_1[[5]], scenario_1[[6]])
start_time <- Sys.time()
SLM_X_1 <- mclapply(prova, SLM_Xreg_Gibbs_parallel,
                  R=500, burnin=1750,
                  dim_cluster=10,
                  mc.cores = 8)
Sys.time() - start_time


start_time <- Sys.time()
SLM_1 <- mclapply(scenario_1, SLM_Gibbs_parallel,
                  R=750, burnin=2000,
                  dim_cluster=10,
                  mc.cores = 8)
Sys.time() - start_time
SLM_2 <- mclapply(scenario_2, SLM_Gibbs_parallel,
                  R=750, burnin=2000,
                  dim_cluster=10,
                  mc.cores = 8)
SLM_3 <- mclapply(scenario_3, SLM_Gibbs_parallel,R=2000, burnin=4000,
                  dim_cluster=10,
                  mc.cores = 6)

start_time <- Sys.time()
SLM_X_1 <- mclapply(scenario_1, SLM_Xreg_Gibbs_parallel,
                  R=500, burnin=1500,
                  dim_cluster=8,
                  mc.cores = 8)
Sys.time() - start_time
SLM_X_2 <- mclapply(scenario_2, SLM_Xreg_Gibbs_parallel,
                    R=500, burnin=1500,
                    dim_cluster=8,
                    mc.cores = 8)
SLM_X_3 <- mclapply(scenario_3, SLM_Xreg_Gibbs_parallel,
                  R=500, burnin=1500,
                  dim_cluster=8,
                  mc.cores = 8)

#############################################################################

# bias and MSE estimation 

bias_funct<-function(data_sim, model_est){
  
  samples= length(model_est)
  
  ATE_Y <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - data_sim[[s]]$simulated_full$Y_0 -
           (model_est[[s]]$Y1_POST_median - model_est[[s]]$Y0_POST_median)))
  ATE_P <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_1 - data_sim[[s]]$simulated_full$P_0 -
           (model_est[[s]]$P1_POST_median - model_est[[s]]$P0_POST_median)))
  bias_P0 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_0 - model_est[[s]]$P0_POST_mean))
  bias_P1 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_1 - model_est[[s]]$P1_POST_mean)) 
  bias_Y0 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_0 - model_est[[s]]$Y0_POST_mean))
  bias_Y1 <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - model_est[[s]]$Y1_POST_mean))
  
  return(list(bias_ATE_Y = ATE_Y, bias_ATE_P = ATE_P, 
              bias_P0 = bias_P0, bias_P1=bias_P1, 
              bias_Y0=bias_Y0, bias_Y1=bias_Y1))
}

bias_funct_SLM<-function(data_sim, model_est){
  
  samples= length(model_est)
  
  ATE_Y <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$Y_1 - data_sim[[s]]$simulated_full$Y_0 -
           (model_est[[s]]$post_Y_1_imp - model_est[[s]]$post_Y_0_imp)))
  ATE_P <- sapply(1:samples, function(s)
    mean(data_sim[[s]]$simulated_full$P_1 - data_sim[[s]]$simulated_full$P_0 -
           (model_est[[s]]$post_P_1_imp - model_est[[s]]$post_P_0_imp)))
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

bias_1 = bias_funct(data_sim = scenario_1, 
                    model_est=model_setting_1)
bias_2_a = bias_funct(data_sim = scenario_2, 
                      model_est=model_setting_2_a)
bias_2_b = bias_funct(data_sim = scenario_2, 
                      model_est=model_setting_2_b)
bias_3_a = bias_funct(data_sim = scenario_3, 
                      model_est=model_setting_3_a)
bias_3_b = bias_funct(data_sim = scenario_3, 
                      model_est=model_setting_3_b)

bias_1_slm = bias_funct_SLM(data_sim = scenario_1, 
                            model_est=SLM_1)
bias_2_slm = bias_funct_SLM(data_sim = scenario_2, 
                            model_est=SLM_2)

bias_1_slm_x = bias_funct_SLM(data_sim = scenario_1, 
                            model_est=SLM_X_1)
bias_2_slm_x = bias_funct_SLM(data_sim = scenario_2, 
                            model_est=SLM_X_2)
bias_3_slm_x = bias_funct_SLM(data_sim = scenario_3, 
                              model_est=SLM_X_3)

sapply(bias_1,function(i) summary(i))
sapply(bias_2_a,function(i) summary(i))
sapply(bias_2_b,function(i) summary(i))
sapply(bias_3_a,function(i) summary(i))
sapply(bias_3_b,function(i) summary(i))


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


#par(mfrow=c(2,2))
#sample_number=5
#plot_var_comparison(sim=scenario_1[[sample_number]]$simulated_full,
#                    estimation=model_setting_1[[sample_number]],
#                    n_scenario=1)


boxplot_allbias(biases=list(bias_1, bias_2_a, bias_2_b, 
                            bias_3_a, bias_3_b))

boxplot_allbias(biases=list(bias_1, bias_1_slm, bias_1_slm_x))
boxplot_allbias(biases=list(bias_2_a,bias_2_b, bias_2_slm, bias_2_slm_x))
boxplot_allbias(biases=list(bias_3_a,bias_3_b, bias_3_slm_x))

#############################################################################
biases=list(bias_1, bias_1_slm_x, 
            bias_2_b,bias_2_slm_x,
            bias_3_b, bias_3_slm_x)
dim_setttings=length(biases)
bias_ATE_Y_all <- cbind(bias=c(sapply(biases, function(j) j$bias_ATE_Y)), 
                        set=rep(1:dim_setttings, each=samples))
bias_ATE_P_all <- cbind(bias=c(sapply(biases, function(j) j$bias_ATE_P)), 
                        set=rep(1:dim_setttings, each=samples))

colour_palette=c("#FAD401","#E09600","#FAD401","#E09600","#FAD401","#E09600")

par(mfrow=c(1,1))
boxplot(bias ~ set, data=bias_ATE_Y_all, 
        main="Y(1)-Y(0)", ylab="bias", xlab="Scenarios",
        ylim=c(-5,5), #500
        col=colour_palette[1:dim_setttings])
abline(h=0, col="#E02700")

boxplot(bias ~ set, data=bias_ATE_P_all, 
        main="P(1)-P(0)", ylab="bias", xlab="Scenarios",
        #ylim=c(-3,3),
        col=colour_palette[1:dim_setttings])
abline(h=0, col="#E02700")

##############################################################################

