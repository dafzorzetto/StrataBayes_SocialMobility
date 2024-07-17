P <- 7      # Num of predictors
n <- 300
c <- 1
samples <- 20

sim_1_z <- function(seed, n){
  
  # set seed for reproducibility
  set.seed(seed)
  cov <- list()
  for(i in 1:P){
    cov[[i]] <- rnorm(n,0,1)
  }
  X <- do.call(cbind, cov)
  X <- cbind(1,X)
  
  reg1 <- ifelse(X[,1]< 0, 1, -1 )
  reg2 <- ifelse(X[,2]< 0, -1, 1 )
  
  Tr <- rbinom(n, 1, pnorm(0.5+reg1+reg2-0.5*abs(X[,3]-1)+1.5*X[,4]*X[,5])) # exposure model
  
  P_err <- rnorm(n, 0, 0.1)
  P_1 <- 0.5*reg1 + 0.5*reg2 + 1*2 + (X[,3]+1) + 1.5*X[,4] - exp(0.3*X[,5]) + (X[,5]) + P_err
  P_0 <- 0.5*reg1 + 0.5*reg2 + (0)*2 + (X[,3]+1) + 1.5*X[,4] - exp(0.3*X[,5]) + 1*(0)*(X[,5]) + P_err
  P_mis <- ifelse(Tr==0, P_1, P_0)
  P_obs <- ifelse(Tr==1, P_1, P_0)
  
  
  S <- P_1-P_0 
  
  Y_err <- rnorm(n, 0, 0.3)
  Y_1 <- 1*reg1+1.5*reg2-(S)*1+2*abs(X[,3]+1)   + 2*X[,4]+ exp(0.5*X[,5])-0.5*abs(X[,6]) - 1*abs(X[,7]+1)+Y_err
  Y_0 <- 1*reg1+1.5*reg2-(S)*(0)+2*abs(X[,3]+1) + 2*X[,4]+ exp(0.5*X[,5])-0.5*abs(X[,6]) - 1*abs(X[,7]+1)+Y_err
  Y_mis <- ifelse(Tr==0, Y_1, Y_0)
  Y_obs <- ifelse(Tr==1, Y_1, Y_0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, 
                Y_obs = Y_obs
    ),
    simulated_full = list(P_0 = P_0,
                          P_1 = P_1,
                          Y_0 = Y_0,
                          Y_1 = Y_1)))
}

sim_2_z <- function(seed, n){
  
  # set seed for reproducibility
  set.seed(seed)
  cov <- list()
  for(i in 1:P){
    cov[[i]] <- rnorm(n,0,1)
  }
  X <- do.call(cbind, cov)
  X <- cbind(1,X)
  
  reg1 <- ifelse(X[,2]< 0, 1, -1 )
  reg2 <- ifelse(X[,3]< 0, -1, 1 )
  
  Tr <- rbinom(n, 1, pnorm(0.5+reg1+reg2-0.5*abs(X[,4]-1)+1.5*X[,5]*X[,6])) # exposure model
  
  P_err <- rnorm(n, 0, 0.1)
  P_1 <- 0.5*reg1 + 0.5*reg2 + 1*2 + (X[,4]+1) + 1.5*X[,5] - exp(0.3*X[,6]) + (X[,6]) + P_err
  P_0 <- 0.5*reg1 + 0.5*reg2 + (0)*2 + (X[,4]+1) + 1.5*X[,5] - exp(0.3*X[,6]) + 1*(0)*(X[,6]) + P_err
  P_mis <- ifelse(Tr==0, P_1, P_0)
  P_obs <- ifelse(Tr==1, P_1, P_0)
  
  # Cluster allocation
  V <- matrix(ifelse(reg1 == 1 & reg2 == 1, 1,
                     ifelse(reg1 == -1 & reg2 == -1, 3, 2)), 
              ncol = 1)
  
  S <- P_1-P_0 
  
  Y_err <- rnorm(n, 0, 0.3)
  Y_1 <- ifelse(V == 1, 4*reg1 + 6*reg2,    
                ifelse(V == 2, 1*reg1 + 1.5*reg2,  
                       (-10)*reg1 + (-15)*reg2)) +  
    -(S)*1 + 2*abs(X[,4] + 1) + 2*X[,5] + exp(0.5*X[,6]) - 0.5*abs(X[,7]) - 1*abs(X[,8] + 1) + Y_err
  Y_0 <- 1*reg1+1.5*reg2-(S)*(0)+2*abs(X[,4]+1) + 2*X[,5]+ exp(0.5*X[,6])-0.5*abs(X[,7]) - 1*abs(X[,8]+1)+Y_err
  Y_mis <- ifelse(Tr==0, Y_1, Y_0)
  Y_obs <- ifelse(Tr==1, Y_1, Y_0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, 
                Y_obs = Y_obs
    ),
    simulated_full = list(P_0 = P_0,
                          P_1 = P_1,
                          Y_0 = Y_0,
                          Y_1 = Y_1)))
}


sim_3_z <- function(seed, n){
  
  # set seed for reproducibility
  set.seed(seed)
  cov <- list()
  for(i in 1:P){
    cov[[i]] <- rnorm(n,0,1)
  }
  X <- do.call(cbind, cov)
  X <- cbind(1,X)
  
  reg1 <- ifelse(X[,2]< 0, 1, -1 )
  reg2 <- ifelse(X[,3]< 0, -1, 1 )
  
  Tr <- rbinom(n, 1, pnorm(0.5+reg1+reg2-0.5*abs(X[,4]-1)+1.5*X[,5]*X[,6])) # exposure model
  
  P_err <- rnorm(n, 0, 0.1)
  P_1 <- 0.5*reg1 + 0.5*reg2 + 1*2 + (X[,4]+1) + 1.5*X[,5] - exp(0.3*X[,6]) + (X[,6]) + P_err
  P_0 <- 0.5*reg1 + 0.5*reg2 + (0)*2 + (X[,4]+1) + 1.5*X[,5] - exp(0.3*X[,6]) + 1*(0)*(X[,6]) + P_err
  P_mis <- ifelse(Tr==0, P_1, P_0)
  P_obs <- ifelse(Tr==1, P_1, P_0)
  
  # Cluster allocation
  V <- matrix(ifelse(reg1 == 1 & reg2 == 1, 1,
                     ifelse(reg1 == -1 & reg2 == -1, 4, 
                            ifelse(reg1 == 1 & reg2 == -1, 2,3))), 
              ncol = 1)
  
  S <- P_1-P_0 
  
  Y_err <- rnorm(n, 0, 0.3)
  Y_1 <- ifelse(V == 1, 10*reg1 + 15*reg2,    
                ifelse(V == 2, 5*reg1 + 0.5*reg2,
                       ifelse(V==3, 10*reg1 + 1*reg2,
                              (10)*reg1 + (15)*reg2))) +  
    -(S)*1 + 2*abs(X[,4] + 1) + 2*X[,5] + exp(0.5*X[,6]) - 0.5*abs(X[,7]) - 1*abs(X[,8] + 1) + Y_err
  Y_0 <- 1*reg1+1.5*reg2-(S)*(0)+2*abs(X[,4]+1) + 2*X[,5]+ exp(0.5*X[,6])-0.5*abs(X[,7]) - 1*abs(X[,8]+1)+Y_err
  Y_mis <- ifelse(Tr==0, Y_1, Y_0)
  Y_obs <- ifelse(Tr==1, Y_1, Y_0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, 
                Y_obs = Y_obs
    ),
    simulated_full = list(P_0 = P_0,
                          P_1 = P_1,
                          Y_0 = Y_0,
                          Y_1 = Y_1)))
}



scenario_1_z <- lapply(1:samples, function(c) sim_1_z(seed = c, n = n))
scenario_2_z <- lapply(1:samples, function(c) sim_2_z(seed = c, n = n))
scenario_3_z <- lapply(1:samples, function(c) sim_3_z(seed = c, n = n))

plot_sim <- function(scenario,title){
  par(mfrow=c(2,2))
  hist(scenario$simulated_full$P_0, breaks = 200, main=NULL, xlab="simulated P0")
  hist(scenario$simulated_full$Y_0, breaks = 200, main=NULL, xlab="simulated Y0")
  hist(scenario$simulated_full$P_1, breaks = 200, main=NULL, xlab="simulated P1")
  hist(scenario$simulated_full$Y_1, breaks = 200, main=NULL, xlab="simulated Y1")
  mtext(title, side=3, outer=TRUE, line=-1, cex=1.5)
}
plot_sim(scenario=scenario_1_z[[1]],title="Scenario 1 Z")
plot_sim(scenario=scenario_2_z[[1]],title="Scenario 2 Z")
plot_sim(scenario=scenario_3_z[[1]],title="Scenario 3 Z")



library(parallel)

Gibbs_parallel_Zigler<- function(data, R=5000, burnin=5000, dim_cluster=10){
  
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


model_1_zigler <- mclapply(scenario_1_z, Gibbs_parallel_Zigler,
                            R=5000, burnin=5000,
                            dim_cluster=10,
                            mc.cores = 5)
model_2_zigler <- mclapply(scenario_2_z, Gibbs_parallel_Zigler,
                           R=5000, burnin=5000,
                           dim_cluster=10,
                           mc.cores = 5)
model_3_zigler <- mclapply(scenario_3_z, Gibbs_parallel_Zigler,
                           R=5000, burnin=5000,
                           dim_cluster=10,
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

bias_1_z = bias_funct(data_sim = scenario_1_z, 
                      model_est=model_1_zigler)
bias_2_z = bias_funct(data_sim = scenario_2_z, 
                      model_est=model_2_zigler)
bias_3_z = bias_funct(data_sim = scenario_3_z, 
                      model_est=model_3_zigler)

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


boxplot_allbias(biases=list(bias_1, bias_2_correct, bias_2_all, bias_1_z,bias_2_z, bias_3_z))


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
cluster_estimated(model_1_zigler)
cluster_estimated(model_2_zigler)
cluster_estimated(model_3_zigler)
