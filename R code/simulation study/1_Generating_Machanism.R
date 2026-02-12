n <- 500 # units for each sample
samples <- 200 # repetition of each setting (number of samples)


### --- scenarios based on first one ---

sim_funct_1 <- function(seed, n) {
  # set seed for reproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectively
  X <- cbind(
    1,
    rbinom(n, 1, 0.4),
    rbinom(n, 1, 0.6),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1)
  )
  
  # treatment
  # logit of a covariates function
  reg_T <- 0 + 0.4 * (X[, 2]) + 0.4 * (X[, 3]) + 0.15 * (X[, 5])
  logit_T <- exp(reg_T) / (1 + exp(reg_T))
  Tr <- rbinom(n, 1, logit_T)
  
  #variabes for post-treatment
  beta_0 <- matrix(c(1, 2, 3, 0.5, 0.1, 0.3), 6)
  beta_1 <- matrix(c(1, 4, 5, 0.5, 0.4, 0.2), 6)
  sigma_0 <- 1
  sigma_1 <- 1
  
  # post-treatment variable
  P_0 <- rnorm(n, X %*% beta_0, sigma_0)
  P_1 <- rnorm(n, X %*% beta_1, sigma_1)
  
  P_obs <- P_0 * (Tr == 0) + P_1 * (Tr == 1)
  P_mis <- P_0 * (Tr == 1) + P_1 * (Tr == 0)
  
  X_0 <- X[Tr == 0, , drop = FALSE]
  X_1 <- X[Tr == 1, , drop = FALSE]
  P_T_0 <- P_0[Tr == 0]
  P_T_1 <- P_1[Tr == 1]
  
  # Cluster allocation
  V <- matrix(ifelse(X[,2] == 1 & X[,3] == 1, 1,
                     ifelse(X[,2] == 0 & X[,3] == 0, 3, 2)), 
              ncol = 1)
  
  # Variables for outcome
  eta_0 <- matrix(c(10,1.5,1.3,2, 2, 0.3, 1.8,
                    1,1.1,0.75,1, 0.2, 0.3, 0.2,
                    -5,0.25,0.1,1, 0.2, 0.1, -1), nrow=3, ncol=7, byrow=TRUE)
  eta_1 <- matrix(c(3,1,1,1,0.3, 0.3, 1.5,
                    0.5,0.5,0.5,0.5,0.3, 0.3, 0.6,
                    -2,0.1,0.1,0.1, 0.4, 0.3, -0.6), nrow=3, ncol=7, byrow=TRUE) #make values smaller
  sigma_y_0 <- matrix(c(1,2,1.5))
  sigma_y_1 <- matrix(c(2,0.5,1))
  
  # Outcome
  #create for every unit a mean vector 
  mean_vector_0 <- rowSums(cbind(X[,1:4],-0.5*abs(X[,5]),exp(0.5*X[,6]),P_0)*eta_0[V,])
  mean_vector_1 <- rowSums(cbind(X[,1:4],-0.5*abs(X[,5]),exp(0.5*X[,6]),(P_1-P_0))*eta_1[V,])
  Y_0 <- rnorm(n,mean_vector_0, sqrt(sigma_y_0[V,]))
  Y_1 <- rnorm(n,mean_vector_1, sqrt(sigma_y_1[V,]))
  
  Y_obs <- Y_0 * (Tr == 0) + Y_1 * (Tr == 1)
  Y_mis <- Y_0 * (Tr == 1) + Y_1 * (Tr == 0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, 
                Y_obs = Y_obs, 
                X_0 = X_0, 
                X_1 = X_1,
                P_T_0 = P_T_0, 
                P_T_1 = P_T_1), # simulated data
    simulated_full = list(P_0 = P_0, P_1 = P_1, Y_0 = Y_0, Y_1 = Y_1),
    parameters = list(
      beta_0 = beta_0, beta_1 = beta_1, # true parameters
      sigma_0 = sigma_0, sigma_1 = sigma_1,
      eta_0 = eta_0, eta_1= eta_1,
      sigma_y_0 = sigma_y_0, sigma_y_1 = sigma_y_1,
      V = V
    )
  ))
}

scenario_1 <- lapply(1:samples, function(c) sim_funct_1(seed = c, n = 500))

### --- scenarios based on second one ---

sim_funct_2 <- function(seed, n) {
  # set seed for reproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectively
  X <- cbind(
    1,
    rbinom(n, 1, 0.4),
    rbinom(n, 1, 0.6),
    rnorm(n, 0, 1),
    rnorm(n, 0, 0.5),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1)
  )
  
  # treatment
  # logit of a covariates function
  reg_T <- 0 + 0.2 * (X[, 2]) + 0.4 * (X[, 3]) + 0.1 * (X[, 7])
  logit_T <- exp(reg_T) / (1 + exp(reg_T))
  Tr <- rbinom(n, 1, logit_T)
  
  #variabes for post-treatment
  beta_0 <- matrix(c(-1, 0.5, 1.5, 0.2, 0.5,
                     0.7, 1 ,-0.5, -1.2), 9)
  beta_1 <- matrix(c(-0.5, 1, 1.8, 0.2, 0.5,
                     0.7, 1.2 ,-0.3, -1), 9)
  sigma_0 <- 0.5
  sigma_1 <- 0.5
  
  # post-treatment variable
  P_0 <- rnorm(n, X[,-c(2:3)] %*% beta_0, sigma_0)
  P_1 <- rnorm(n, X[,-c(2:3)] %*% beta_1, sigma_1)
  
  P_obs <- P_0 * (Tr == 0) + P_1 * (Tr == 1)
  P_mis <- P_0 * (Tr == 1) + P_1 * (Tr == 0)
  
  X_0 <- X[Tr == 0, , drop = FALSE]
  X_1 <- X[Tr == 1, , drop = FALSE]
  P_T_0 <- P_0[Tr == 0]
  P_T_1 <- P_1[Tr == 1]
  
  # Cluster allocation
  V <- matrix(ifelse(X[,2] == 1 & X[,3] == 1, 1,
                     ifelse(X[,2] == 0 & X[,3] == 0, 3, 2)), 
              ncol = 1)
  
  # Variables for outcome
  eta_0 <- matrix(c(10,1.5,1.3,0.1, 0.4,0.1, 0.2,-0.4,0.3,2,
                    1,1.1,0.75,0.1,0.1, 0.2, 0.2,-0.4,0.3,1,
                    -5,0.25,0.1,0.5,0.1, 0.2, 0.4,-0.4,-0.3,-1), nrow=3, ncol=10, byrow=TRUE)
  eta_1 <- matrix(c(10,1,1,0.6, 0.1,0.1, 0.2,-0.4,0.3,2.5,
                    0,0.8,0.5,0.5,0.1, 0.2, 0.4,-0.4,0.3,0.5,
                    -5,0.5,0.25,0.2,0.1, 0.2, 0.7,-0.4,0.4,-2), nrow=3, ncol=10, byrow=TRUE) #make values smaller
  sigma_y_0 <- matrix(c(0.5,0.5,0.5))
  sigma_y_1 <- matrix(c(0.5,0.5,0.5))
  
  # Outcome
  #create for every unit a mean vector 
  mean_vector_0 <- rowSums(cbind(X[,c(1,4,5,6,7)],abs(X[,8]+2),exp(X[,9]*0.2),-0.1*abs(X[,10]),exp(0.1*X[,11]),P_0)*eta_0[V,])
  mean_vector_1 <- rowSums(cbind(X[,c(1,4,5,6,7)],abs(X[,8]+1.5),exp(X[,9]*0.1),-0.3*abs(X[,10]),exp(0.2*X[,11]),(P_1-P_0))*eta_1[V,])
  Y_0 <- rnorm(n,mean_vector_0, sqrt(sigma_y_0[V,]))
  Y_1 <- rnorm(n,mean_vector_1, sqrt(sigma_y_1[V,]))
  
  Y_obs <- Y_0 * (Tr == 0) + Y_1 * (Tr == 1)
  Y_mis <- Y_0 * (Tr == 1) + Y_1 * (Tr == 0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, 
                Y_obs = Y_obs, 
                X_0 = X_0, 
                X_1 = X_1,
                P_T_0 = P_T_0, 
                P_T_1 = P_T_1), # simulated data
    simulated_full = list(P_0 = P_0, P_1 = P_1, Y_0 = Y_0, Y_1 = Y_1),
    parameters = list(
      beta_0 = beta_0, beta_1 = beta_1, # true parameters
      sigma_0 = sigma_0, sigma_1 = sigma_1,
      eta_0 = eta_0, eta_1= eta_1,
      sigma_y_0 = sigma_y_0, sigma_y_1 = sigma_y_1,
      V = V
    )
  ))
}

scenario_2 <- lapply(1:samples, function(c) sim_funct_2(seed = c, n = 500))

sim_funct_3 <- function(seed, n) {
  # set seed for reproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectively
  X <- cbind(
    1,
    rbinom(n, 1, 0.4),
    rbinom(n, 1, 0.6),
    rnorm(n, 0, 1),
    rnorm(n, 0, 0.5),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 1),
    rnorm(n, 0, 0.2),
    rnorm(n, 0, 0.25),
    rnorm(n, 0, 0.50),
    rnorm(n, 0, 0.75)
  )
  
  # treatment
  # logit of a covariates function
  reg_T <- 0 + 0.2 * (X[, 2]) + 0.4 * (X[, 3]) + 0.1 * (X[, 7])
  logit_T <- exp(reg_T) / (1 + exp(reg_T))
  Tr <- rbinom(n, 1, logit_T)
  
  #variabes for post-treatment
  beta_0 <- matrix(c(-1, 0.5, 1.5, 0.2, 0.5,
                     0.7, 1 ,-0.5, -1.2, 0.1,
                     0.1, 0.1, 0.1), 13)
  beta_1 <- matrix(c(-0.5, 1, 1.8, 0.2, 0.5,
                     0.7, 1.2 ,-0.3, -1, 0.1,
                     0.1, 0.1, 0.1), 13)
  sigma_0 <- 0.5
  sigma_1 <- 0.5
  
  # post-treatment variable
  P_0 <- rnorm(n, X[,-c(2:3)] %*% beta_0, sigma_0)
  P_1 <- rnorm(n, X[,-c(2:3)] %*% beta_1, sigma_1)
  
  P_obs <- P_0 * (Tr == 0) + P_1 * (Tr == 1)
  P_mis <- P_0 * (Tr == 1) + P_1 * (Tr == 0)
  
  X_0 <- X[Tr == 0, , drop = FALSE]
  X_1 <- X[Tr == 1, , drop = FALSE]
  P_T_0 <- P_0[Tr == 0]
  P_T_1 <- P_1[Tr == 1]
  
  # Cluster allocation
  V <- matrix(ifelse(X[,2] == 1 & X[,3] == 1, 1,
                     ifelse(X[,2] == 0 & X[,3] == 0, 3, 2)), 
              ncol = 1)
  
  # Variables for outcome
  eta_0 <- matrix(c(10,1.5,1.3,0.1,0.4,0.1,0.2,-0.4,0.3,0.1,0.1,0.1,0.1,2,
                    1,1.1,0.75,0.1,0.1, 0.2, 0.2,-0.4,0.3,0.1,0.1,0.1,0.1,1,
                    -5,0.25,0.1,0.5,0.1, 0.2, 0.4,-0.4,-0.3,0.1,0.1,0.1,0.1,-1), nrow=3, ncol=14, byrow=TRUE)
  eta_1 <- matrix(c(10,1,1,0.6, 0.1,0.1, 0.2,-0.4,0.3,0.1,0.1,0.1,0.1,2.5,
                    0,0.8,0.5,0.5,0.1, 0.2, 0.4,-0.4,0.3,0.1,0.1,0.1,0.1,0.5,
                    -5,0.5,0.25,0.2,0.1, 0.2, 0.7,-0.4,0.4,0.1,0.1,0.1,0.1,-2), nrow=3, ncol=14, byrow=TRUE)
  sigma_y_0 <- matrix(c(0.5,0.5,0.5))
  sigma_y_1 <- matrix(c(0.5,0.5,0.5))
  
  # Outcome
  #create for every unit a mean vector 
  mean_vector_0 <- rowSums(cbind(X[,c(1,4,5,6,7)],abs(X[,8]+2),exp(X[,9]*0.2),
                                 -0.1*abs(X[,10]),exp(0.1*X[,11]),X[,c(12,13,14,15)],
                                 P_0)*eta_0[V,])
  mean_vector_1 <- rowSums(cbind(X[,c(1,4,5,6,7)],abs(X[,8]+1.5),exp(X[,9]*0.1),
                                 -0.3*abs(X[,10]),exp(0.2*X[,11]),X[,c(12,13,14,15)],
                                 (P_1-P_0))*eta_1[V,])
  Y_0 <- rnorm(n,mean_vector_0, sqrt(sigma_y_0[V,]))
  Y_1 <- rnorm(n,mean_vector_1, sqrt(sigma_y_1[V,]))
  
  Y_obs <- Y_0 * (Tr == 0) + Y_1 * (Tr == 1)
  Y_mis <- Y_0 * (Tr == 1) + Y_1 * (Tr == 0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, 
                Y_obs = Y_obs, 
                X_0 = X_0, 
                X_1 = X_1,
                P_T_0 = P_T_0, 
                P_T_1 = P_T_1), # simulated data
    simulated_full = list(P_0 = P_0, P_1 = P_1, Y_0 = Y_0, Y_1 = Y_1),
    parameters = list(
      beta_0 = beta_0, beta_1 = beta_1, # true parameters
      sigma_0 = sigma_0, sigma_1 = sigma_1,
      eta_0 = eta_0, eta_1= eta_1,
      sigma_y_0 = sigma_y_0, sigma_y_1 = sigma_y_1,
      V = V
    )
  ))
}

scenario_3 <- lapply(1:samples, function(c) sim_funct_3(seed = c, n = 300))


# plots ####
plot_potential_P<-function(sim_P0, sim_P1, sim_V, n_scenario){
  plot(sim_P0,sim_P1, pch=19, main=paste0("scenario ", n_scenario))
  points(sim_P0[sim_V==2],
         sim_P1[sim_V==2], pch=19,
         col=2)
  points(sim_P0[sim_V==3],
         sim_P1[sim_V==3], pch=19,
         col=3)
}
plot_potential_Y<-function(sim_Y0, sim_Y1, sim_V, n_scenario){
  plot(sim_Y0,sim_Y1, pch=19, main=paste0("scenario ", n_scenario))
  points(sim_Y0[sim_V==2],
         sim_Y1[sim_V==2], pch=19,
         col=2)
  points(sim_Y0[sim_V==3],
         sim_Y1[sim_V==3], pch=19,
         col=3)
}

### scenario 1
par(mfrow=c(2,3))
plot_potential_P(scenario_1_3_2_X[[1]]$simulated_full$P_0,
                 scenario_1_3_2_X[[1]]$simulated_full$P_1,
                 scenario_1_3_2_X[[1]]$parameters$V,
                 n_scenario="1.3.2x")
###scenario 2
plot_potential_P(scenario_2_3_2_X[[1]]$simulated_full$P_0,
                 scenario_2_3_2_X[[1]]$simulated_full$P_1,
                 scenario_2_3_2_X[[1]]$parameters$V,
                 n_scenario="2.3.2x")
plot_potential_P(scenario_2_3_3_X[[1]]$simulated_full$P_0,
                 scenario_2_3_3_X[[1]]$simulated_full$P_1,
                 scenario_2_3_3_X[[1]]$parameters$V,
                 n_scenario="2.3.3x")
### scenario 1
plot_potential_Y(scenario_1_3_2_X[[1]]$simulated_full$Y_0,
                 scenario_1_3_2_X[[1]]$simulated_full$Y_1,
                 scenario_1_3_2_X[[1]]$parameters$V,
                 n_scenario="1.3.2x")
### scenario 2
plot_potential_Y(scenario_2_3_2_X[[1]]$simulated_full$Y_0,
                 scenario_2_3_2_X[[1]]$simulated_full$Y_1,
                 scenario_2_3_2_X[[1]]$parameters$V,
                 n_scenario="2.3.2x")
plot_potential_Y(scenario_2_3_3_X[[1]]$simulated_full$Y_0,
                 scenario_2_3_3_X[[1]]$simulated_full$Y_1,
                 scenario_2_3_3_X[[1]]$parameters$V,
                 n_scenario="2.3.3x")

### histogram simulation function ###
hist_simulation <- function(scenario, title){
  par(mfrow=c(2,2))
  hist(scenario[[1]]$simulated_full$P_0, breaks = 200, main="P0", xlab="simulated P0")
  hist(scenario[[1]]$simulated_full$Y_0, breaks = 200, main="Y0", xlab="simulated Y0")
  hist(scenario[[1]]$simulated_full$P_1, breaks = 200, main="P1", xlab="simulated P1")
  hist(scenario[[1]]$simulated_full$Y_1, breaks = 200, main="Y1", xlab="simulated Y1")
  mtext(title, side=3, line=-2, outer=TRUE, cex=1.5)
}


hist_simulation(scenario_1_3_2_X, "Scenario 1.3.2 X")

hist_simulation(scenario_2_3_2_X, "Scenario 2.3.2 X")

hist_simulation(scenario_2_3_3_X, "Scenario 2.3.3 X")


#hist_simulation_condTR(scenario_1)
#hist_simulation_condTR(scenario_1_3_X)
#hist_simulation_condTR(scenario_1_2_X)

#hist_simulation_condTR(scenario_2)
#hist_simulation_condTR(scenario_2_3_X)
#hist_simulation_condTR(scenario_2_2_X)




