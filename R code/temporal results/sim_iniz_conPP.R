n <- 500 # units for each sample
samples <- 20 # repetition of each setting (number of samples)


### --- scenarios based on first one ---
sim_bin_2 <- function(seed, n) {
  # set seed for reproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectively
  X <- cbind(
    1,
    rbinom(n, 1, 0.4),
    rbinom(n, 1, 0.6)
  )
  
  # treatment
  # logit of a covariates function
  reg_T <- 0 + 0.4 * (X[, 2]) + 0.6 * (X[, 3])
  logit_T <- exp(reg_T) / (1 + exp(reg_T))
  Tr <- rbinom(n, 1, logit_T)
  
  #variabes for post-treatment
  beta_0 <- matrix(c(1, 2, 3), 3)
  beta_1 <- matrix(c(1, 4, 5), 3)
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
  eta_0 <- matrix(c(10,1.5,1.3,2,1,1.1,0.75,1,-5,0.25,0.1,1), nrow=3, ncol=4, byrow=TRUE)
  eta_1 <- matrix(c(1,1,1,1,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1), nrow=3, ncol=4, byrow=TRUE) #make values smaller
  sigma_y_0 <- matrix(c(1,2,1.5))
  sigma_y_1 <- matrix(c(2,0.5,1))
  
  # Outcome
  #create for every unit a mean vector 
  mean_vector_0 <- rowSums(cbind(X,P_0)*eta_0[V,])
  mean_vector_1 <- rowSums(cbind(X,(P_1-P_0)^2)*eta_1[V,])
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

sim_bin_3 <- function(seed, n) {
  # set seed for reproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectively
  X <- cbind(
    1,
    rbinom(n, 1, 0.4),
    rbinom(n, 1, 0.6)
  )
  
  # treatment
  # logit of a covariates function
  reg_T <- 0 + 0.4 * (X[, 2]) + 0.6 * (X[, 3])
  logit_T <- exp(reg_T) / (1 + exp(reg_T))
  Tr <- rbinom(n, 1, logit_T)
  
  #variabes for post-treatment
  beta_0 <- matrix(c(1, 2, 3), 3)
  beta_1 <- matrix(c(1, 4, 5), 3)
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
  eta_0 <- matrix(c(10,1.5,1.3,2,1,1.1,0.75,1,-5,0.25,0.1,1), nrow=3, ncol=4, byrow=TRUE)
  eta_1 <- matrix(c(1,1,1,1,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1), nrow=3, ncol=4, byrow=TRUE) #make values smaller
  sigma_y_0 <- matrix(c(1,2,1.5))
  sigma_y_1 <- matrix(c(2,0.5,1))
  
  # Outcome
  #create for every unit a mean vector 
  mean_vector_0 <- rowSums(cbind(X,P_0)*eta_0[V,])
  mean_vector_1 <- rowSums(cbind(X,(P_1-P_0))*eta_1[V,])
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

scenario_1_2 <- lapply(1:samples, function(c) sim_bin_2(seed = c, n = n))
scenario_1_3 <- lapply(1:samples, function(c) sim_bin_3(seed = c, n = n))



### --- scenarios based on second one ---
sim_bin_cont_2 <- function(seed, n) {
  # set seed for reproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectively
  X <- cbind(
    1,
    rbinom(n, 1, 0.4),
    rbinom(n, 1, 0.6),
    rnorm(n, 0, 1),
    rnorm(n, 0, 0.5)
  )
  
  # treatment
  # logit of a covariates function
  reg_T <- 0 + 0.4 * (X[, 2]) + 0.6 * (X[, 3])
  logit_T <- exp(reg_T) / (1 + exp(reg_T))
  Tr <- rbinom(n, 1, logit_T)
  
  #variabes for post-treatment
  beta_0 <- matrix(c(-1, 0.5, 1.5), 3)
  beta_1 <- matrix(c(-0.5, 1, 1.8), 3)
  sigma_0 <- 0.5
  sigma_1 <- 0.5
  
  # post-treatment variable
  P_0 <- rnorm(n, X[,c(1,4,5)] %*% beta_0, sigma_0)
  P_1 <- rnorm(n, X[,c(1,4,5)] %*% beta_1, sigma_1)
  
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
  eta_0 <- matrix(c(10,1.5,1.3,2,1,1.1,0.75,1,-5,0.25,0.1,1), nrow=3, ncol=4, byrow=TRUE)
  eta_1 <- matrix(c(10,1,1,1,0,0.8,0.5,0.5,-5,0.5,0.25,0.25), nrow=3, ncol=4, byrow=TRUE) #make values smaller
  sigma_y_0 <- matrix(c(0.5,0.5,0.5))
  sigma_y_1 <- matrix(c(0.5,0.5,0.5))
  
  # Outcome
  #create for every unit a mean vector 
  mean_vector_0 <- rowSums(cbind(X[,c(1,4,5)],P_0)*eta_0[V,])
  mean_vector_1 <- rowSums(cbind(X[,c(1,4,5)],(P_1-P_0)^2)*eta_1[V,])
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

sim_bin_cont_3 <- function(seed, n) {
  # set seed for reproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectively
  X <- cbind(
    1,
    rbinom(n, 1, 0.4),
    rbinom(n, 1, 0.6),
    rnorm(n, 0, 1),
    rnorm(n, 0, 0.5)
  )
  
  # treatment
  # logit of a covariates function
  reg_T <- 0 + 0.4 * (X[, 2]) + 0.6 * (X[, 3])
  logit_T <- exp(reg_T) / (1 + exp(reg_T))
  Tr <- rbinom(n, 1, logit_T)
  
  #variabes for post-treatment
  beta_0 <- matrix(c(-1, 0.5, 1.5), 3)
  beta_1 <- matrix(c(-0.5, 1, 1.8), 3)
  sigma_0 <- 0.5
  sigma_1 <- 0.5
  
  # post-treatment variable
  P_0 <- rnorm(n, X[,c(1,4,5)] %*% beta_0, sigma_0)
  P_1 <- rnorm(n, X[,c(1,4,5)] %*% beta_1, sigma_1)
  
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
  eta_0 <- matrix(c(10,1.5,1.3,2,1,1.1,0.75,1,-5,0.25,0.1,1), nrow=3, ncol=4, byrow=TRUE)
  eta_1 <- matrix(c(10,1,1,1,0,0.8,0.5,0.5,-5,0.5,0.25,0.25), nrow=3, ncol=4, byrow=TRUE) #make values smaller
  sigma_y_0 <- matrix(c(0.5,0.5,0.5))
  sigma_y_1 <- matrix(c(0.5,0.5,0.5))
  
  # Outcome
  #create for every unit a mean vector 
  mean_vector_0 <- rowSums(cbind(X[,c(1,4,5)],P_0)*eta_0[V,])
  mean_vector_1 <- rowSums(cbind(X[,c(1,4,5)],(P_1-P_0))*eta_1[V,])
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

scenario_2_2 <- lapply(1:samples, function(c) sim_bin_cont_2(seed = c, n = n))
scenario_2_3 <- lapply(1:samples, function(c) sim_bin_cont_3(seed = c, n = n))


############################################################################

par(mfrow=c(2,3))
plot_potential_P(scenario_1[[1]]$simulated_full$P_0,
                 scenario_1[[1]]$simulated_full$P_1,
                 scenario_1[[1]]$parameters$V,
                 n_scenario=1)
plot_potential_P(scenario_1_3[[1]]$simulated_full$P_0,
                 scenario_1_3[[1]]$simulated_full$P_1,
                 scenario_1_3[[1]]$parameters$V,
                 n_scenario=1.3)
plot_potential_P(scenario_1_2[[1]]$simulated_full$P_0,
                 scenario_1_2[[1]]$simulated_full$P_1,
                 scenario_1_2[[1]]$parameters$V,
                 n_scenario=1.2)
plot_potential_Y(scenario_1[[1]]$simulated_full$Y_0,
                 scenario_1[[1]]$simulated_full$Y_1,
                 scenario_1[[1]]$parameters$V,
                 n_scenario=1)
plot_potential_Y(scenario_1_3[[1]]$simulated_full$Y_0,
                 scenario_1_3[[1]]$simulated_full$Y_1,
                 scenario_1_3[[1]]$parameters$V,
                 n_scenario=1.3)
plot_potential_Y(scenario_1_2[[1]]$simulated_full$Y_0,
                 scenario_1_2[[1]]$simulated_full$Y_1,
                 scenario_1_2[[1]]$parameters$V,
                 n_scenario=1.2)


plot_potential_P(scenario_2[[1]]$simulated_full$P_0,
                 scenario_2[[1]]$simulated_full$P_1,
                 scenario_2[[1]]$parameters$V,
                 n_scenario=2)
plot_potential_P(scenario_2_3[[1]]$simulated_full$P_0,
                 scenario_2_3[[1]]$simulated_full$P_1,
                 scenario_2_3[[1]]$parameters$V,
                 n_scenario=2.3)
plot_potential_P(scenario_2_2[[1]]$simulated_full$P_0,
                 scenario_2_2[[1]]$simulated_full$P_1,
                 scenario_2_2[[1]]$parameters$V,
                 n_scenario=2.2)
plot_potential_Y(scenario_2[[1]]$simulated_full$Y_0,
                 scenario_2[[1]]$simulated_full$Y_1,
                 scenario_2[[1]]$parameters$V,
                 n_scenario=2)
plot_potential_Y(scenario_2_3[[1]]$simulated_full$Y_0,
                 scenario_2_3[[1]]$simulated_full$Y_1,
                 scenario_2_3[[1]]$parameters$V,
                 n_scenario=2.3)
plot_potential_Y(scenario_2_2[[1]]$simulated_full$Y_0,
                 scenario_2_2[[1]]$simulated_full$Y_1,
                 scenario_2_2[[1]]$parameters$V,
                 n_scenario=2.2)


hist_simulation(scenario_1)
hist_simulation(scenario_1_3)
hist_simulation(scenario_1_2)

hist_simulation(scenario_2)
hist_simulation(scenario_2_3)
hist_simulation(scenario_2_2)

hist_simulation_condTR(scenario_1)
hist_simulation_condTR(scenario_1_3)
hist_simulation_condTR(scenario_1_2)

hist_simulation_condTR(scenario_2)
hist_simulation_condTR(scenario_2_3)
hist_simulation_condTR(scenario_2_2)

