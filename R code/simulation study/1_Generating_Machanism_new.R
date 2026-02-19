############################################################################
#       --- simulation dataset ---
############################################################################

n <- 500             # units for each sample
samples <- 100       # repetition of each setting (number of samples)


### --- functions for the simulations ---

# Principal stratification : scenarios 1 and 2
# Mediation : scenarios 3 and 4

Cov_tr <- function(n){
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectively
  X <- cbind(rbinom(n, 1, 0.4),
             rbinom(n, 1, 0.6),
             rnorm(n, 0, 1),
             rnorm(n, 0, 1),
             rnorm(n, 0, 1))
  
  # treatment
  # logit of a covariates function
  reg_T <- -0.1 + 0.3 * (X[, 2]) + 0.2 * (X[, 3]) + 0.1 * (X[, 5])
  logit_T <- exp(reg_T) / (1 + exp(reg_T))
  Tr <- rbinom(n, 1, logit_T)
  
  return(list(X = X, Tr = Tr))
}

sim_funct_1 <- function(seed, n) {
  
  # set seed for reproducibility
  set.seed(seed)
  
  # confounders and tretment
  sim_cov_tr <- Cov_tr(n)
  X <- sim_cov_tr$X
  Tr <- sim_cov_tr$Tr
  
  # post-treatment variable
  X_non_linear <- cbind(1, X[,1], X[,1]*X[,2], sin(X[,3]), X[,4], exp(X[,5]-1))
  alpha <- c(1, 2, -3, 2, -1, 0.5)
  P_0 <- rnorm(n, X_non_linear %*% alpha, 1)
  P_1 <- rnorm(n, cbind(X_non_linear, Tr) %*% c(alpha, 2), 1.3)
  
  P_obs <- P_0 * (Tr == 0) + P_1 * (Tr == 1)
  P_mis <- P_0 * (Tr == 1) + P_1 * (Tr == 0)
  
  # Outcome
  #create for every unit a mean vector 
  mean_Y_com <- cbind(1, X[,1:3],abs(X[,4]),exp(0.5*X[,5]))
  gamma <- c(-2, 1, -1, 2, -1, 0.8)
  Y_0 <- rnorm(n, cbind(mean_Y_com, P_0) %*% c(gamma, 2), 1)
  Y_1 <- rnorm(n, cbind(mean_Y_com, Tr, P_0, P_1) %*% c(gamma, 1.5, 2, 0.5), 1)
  
  Y_obs <- Y_0 * (Tr == 0) + Y_1 * (Tr == 1)
  Y_mis <- Y_0 * (Tr == 1) + Y_1 * (Tr == 0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, Y_obs = Y_obs),
    simulated_full = list(P_0 = P_0, P_1 = P_1, Y_0 = Y_0, Y_1 = Y_1)
  ))
}

sim_funct_2 <- function(seed, n) {
  # set seed for reproducibility
  set.seed(seed)
  
  # confounders and tretment
  sim_cov_tr <- Cov_tr(n)
  X <- sim_cov_tr$X
  Tr <- sim_cov_tr$Tr
  
  # post-treatment variable
  X_non_linear <- cbind(1, X[,1], X[,1]*(1-X[,2]), sin(2*X[,3]), abs(X[,4]), exp(X[,5]-2))
  alpha <- c(-2, 1.5, -3, 1, -1.5, 2)
  P_0 <- rnorm(n, X_non_linear %*% alpha, 1)
  P_1 <- rnorm(n, cbind(X_non_linear, Tr) %*% c(alpha, 1), 1)
  
  P_obs <- P_0 * (Tr == 0) + P_1 * (Tr == 1)
  P_mis <- P_0 * (Tr == 1) + P_1 * (Tr == 0)
  
  # Outcome
  #create for every unit a mean vector 
  mean_Y_com <- cbind(1, X[,1:2], X[,3]*X[,4], abs(X[,4]*0.5),exp(X[,5]))
  gamma <- c(-1, 1.5, -1, 0.5, -1, 0.3)
  Y_0 <- rnorm(n, cbind(mean_Y_com) %*% gamma, 1)
  Y_1 <- rnorm(n, cbind(mean_Y_com, Tr, P_1 - P_0) %*% c(gamma, 1, 1), 1)
  
  Y_obs <- Y_0 * (Tr == 0) + Y_1 * (Tr == 1)
  Y_mis <- Y_0 * (Tr == 1) + Y_1 * (Tr == 0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, Y_obs = Y_obs),
    simulated_full = list(P_0 = P_0, P_1 = P_1, Y_0 = Y_0, Y_1 = Y_1)
  ))
}

sim_funct_3 <- function(seed, n) {
  
  # set seed for reproducibility
  set.seed(seed)
  
  # confounders and tretment
  sim_cov_tr <- Cov_tr(n)
  X <- sim_cov_tr$X
  Tr <- sim_cov_tr$Tr
  
  # post-treatment variable
  X_non_linear <- cbind(1, X[,1:2], X[,1]*X[,4], sin(X[,5]))
  alpha <- c(-2, 1.5, -3, 1, -1.5)
  P_0 <- rnorm(n, X_non_linear %*% alpha, 1)
  P_1 <- rnorm(n, cbind(X_non_linear, Tr) %*% c(alpha, 1), 1)
  
  P_obs <- P_0 * (Tr == 0) + P_1 * (Tr == 1)
  P_mis <- P_0 * (Tr == 1) + P_1 * (Tr == 0)
  
  # Outcome
  #create for every unit a mean vector 
  mean_Y_com <- cbind(1, X[,1:3], abs(X[,4]*0.5),exp(X[,5]))
  gamma <- c(-1, 1.5, -1, 0.5, -1, 0.3)
  Y_0 <- rnorm(n, cbind(mean_Y_com, P_0) %*% c(gamma,1), 1)
  Y_1 <- rnorm(n, cbind(mean_Y_com, Tr, P_1) %*% c(gamma, 1, 0.7), 1)
  
  Y_obs <- Y_0 * (Tr == 0) + Y_1 * (Tr == 1)
  Y_mis <- Y_0 * (Tr == 1) + Y_1 * (Tr == 0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, Y_obs = Y_obs),
    simulated_full = list(P_0 = P_0, P_1 = P_1, Y_0 = Y_0, Y_1 = Y_1)
  ))
}

sim_funct_4 <- function(seed, n) {
  
  # set seed for reproducibility
  set.seed(seed)
  
  # confounders and tretment
  sim_cov_tr <- Cov_tr(n)
  X <- sim_cov_tr$X
  Tr <- sim_cov_tr$Tr
  
  # post-treatment variable
  X_non_linear <- cbind(1, X[,1:2], X[,3]*X[,4], sin(X[,5]))
  alpha <- c(-2, 1.5, -3, 1, -1.5)
  P_0 <- rnorm(n, X_non_linear %*% alpha, 1)
  P_1 <- rnorm(n, cbind(X_non_linear, Tr) %*% c(alpha, 1), 1)
  
  P_obs <- P_0 * (Tr == 0) + P_1 * (Tr == 1)
  P_mis <- P_0 * (Tr == 1) + P_1 * (Tr == 0)
  
  # Outcome
  #create for every unit a mean vector 
  mean_Y_com <- cbind(1, X[,1:3], sin(X[,4]*0.5),exp(X[,5]))
  gamma <- c(-1, 1.5, -1, 0.5, -1, 0.3)
  Y_0 <- rnorm(n, cbind(mean_Y_com, P_0, P_0*X[,1]) %*% c(gamma,1, 0.2), 1)
  Y_1 <- rnorm(n, cbind(mean_Y_com, Tr, P_1, P_1*X[,2]) %*% c(gamma, 1, 1, 0.2), 1)
  
  Y_obs <- Y_0 * (Tr == 0) + Y_1 * (Tr == 1)
  Y_mis <- Y_0 * (Tr == 1) + Y_1 * (Tr == 0)
  
  # save all the information as a list
  return(list(
    data = list(X = X, Tr = Tr, 
                P_obs = P_obs, Y_obs = Y_obs),
    simulated_full = list(P_0 = P_0, P_1 = P_1, Y_0 = Y_0, Y_1 = Y_1)
  ))
}

scenario_1 <- lapply(1:samples, function(c) sim_funct_1(seed = c, n = 500))
scenario_2 <- lapply(1:samples, function(c) sim_funct_2(seed = c, n = 500))
scenario_3 <- lapply(1:samples, function(c) sim_funct_3(seed = c, n = 500))
scenario_4 <- lapply(1:samples, function(c) sim_funct_4(seed = c, n = 300))

save(scenario_1, file = "scenario_1.RData")
save(scenario_2, file = "scenario_2.RData")
save(scenario_3, file = "scenario_3.RData")
save(scenario_4, file = "scenario_4.RData")

### histogram simulation function ###
hist_simulation <- function(scenario, title){
  par(mfrow=c(2,2))
  hist(scenario[[1]]$simulated_full$P_0, breaks = 200, main="P0", xlab="simulated P0")
  hist(scenario[[1]]$simulated_full$Y_0, breaks = 200, main="Y0", xlab="simulated Y0")
  hist(scenario[[1]]$simulated_full$P_1, breaks = 200, main="P1", xlab="simulated P1")
  hist(scenario[[1]]$simulated_full$Y_1, breaks = 200, main="Y1", xlab="simulated Y1")
  mtext(title, side=3, line=-2, outer=TRUE, cex=1.5)
}


hist_simulation(scenario_1, "Scenario 1")
hist_simulation(scenario_2, "Scenario 2")
hist_simulation(scenario_3, "Scenario 3")
hist_simulation(scenario_4, "Scenario 4")


