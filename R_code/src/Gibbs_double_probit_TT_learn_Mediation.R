library(mvtnorm)
library(matrixStats)
library(truncnorm)
library(gtools)         # Dirichlet 

Double_probit_TT_learn_Mediation <- function(Outcome, Interm, Treat, Covariates, 
                                            R_tot, R_burnin, L_max, seed=111, 
                                            method = "population"){
  
  set.seed(seed)
  sample_population = 1000
  
  ####################################################################
  #   ---    pre-process simulated data    ---
  
  n_units =length(Treat)
  if (length(Outcome)!=n_units){
    return(print("check dimension of outcome or treatment vector"))
  }
  if (length(Interm)!=n_units){
    return(print("check dimension of intermediate variable or treatment vector"))
  }
  Covariates <- as.matrix(Covariates)
  if (dim(Covariates)[1]!=n_units){
    return(print("check dimension of covariate matrix"))
  }
  
  # regression matrices
  which_Tr_0 = which(Treat==0)
  which_Tr_1 = which(Treat==1)
  n_units_Tr1 = length(which_Tr_1)
  n_units_Tr0 = length(which_Tr_0)
  
  Interm_0 = Interm[which_Tr_0]
  Interm_1 = Interm[which_Tr_1]
  
  interm_Cov = cbind(1, Covariates)
  interm_Cov_0 = cbind(1, Covariates[which_Tr_0,])
  interm_Cov_1 = cbind(1, Covariates[which_Tr_1,])
  
  weight_Cov = cbind(1, Covariates)
  weight_Cov_0 = weight_Cov[which_Tr_0,]
  weight_Cov_1 = weight_Cov[which_Tr_1,]
  
  Outcome_0 = Outcome[which_Tr_0]
  Outcome_1 = Outcome[which_Tr_1]
  
  outcome_Cov = cbind(interm_Cov, Interm) 
  outcome_Cov_0 = outcome_Cov[which_Tr_0,]
  outcome_Cov_1 = outcome_Cov[which_Tr_1,]
  
  dim_W_cov = ncol(weight_Cov)
  dim_I.tr_cov = ncol(interm_Cov_0)
  dim_O.tr_cov = ncol(outcome_Cov)
  
  Cov_Imp_outcome_Tr0_Interm0 <- Cov_Imp_outcome_Tr0_Interm1 <- outcome_Cov
  Cov_Imp_outcome_Tr1_Interm0 <- Cov_Imp_outcome_Tr1_Interm1 <- outcome_Cov
  
  ####################################################################
  #   ---    functions    ---
  
  stick_breaking_wights <- function(Cov, eta){
    Pnorm_regression = pnorm(Cov %*% eta)
    Pnorm_regression[, L_max] = 1
    log_Pnorm = log(Pnorm_regression)  # logarithm better control the zeros weights
    log_1_P = log(1-Pnorm_regression)[, - L_max]
    cumsum_log1_P = cbind(0, t(apply(log_1_P, 1, cumsum)))
    
    weights = exp(log_Pnorm + cumsum_log1_P)
    return(weights)
  }
  
  Z_sampling <- function(weights, allocation){
    
    denom_regression = cbind(1,1/(1-t(apply(weights[,-L_max], 1, cumsum))))
    phi_regression = weights/denom_regression
    phi_regression[which(phi_regression==1)] <- 1 - .Machine$double.eps
    phi_regression[which(phi_regression<.Machine$double.eps)] <-  .Machine$double.eps
    mean_regression = qnorm(phi_regression)
    
    n_ = length(allocation)
    lim_inf <- matrix(-Inf, nrow=n_, ncol=L_max)
    lim_sup <- matrix(0, nrow=n_, ncol=L_max)
    for(i in 1:n_){
      all_i = allocation[i]
      seq_sup <- if(all_i == L_max) L_max else (all_i+1):L_max
      lim_inf[i, seq_sup] <- NA
      lim_inf[i, all_i] <- 0
      lim_sup[i, seq_sup] <- NA
      lim_sup[i, all_i] <- Inf
    }
    Z_temp = sapply(1:n_, function(i) rtruncnorm(L_max, a = lim_inf[i,], b = lim_sup[i,], mean = mean_regression[i,], sd = 1))
    return(Z_temp)
  }
  
  weights_estimation_post <- function(omega, matrix_reg, par_reg, sigma_reg, out_vec){
    
    n_units_ = nrow(matrix_reg)
    matrix_regression_all = matrix_reg%*%par_reg
    log_gauss = do.call(rbind, lapply(1:n_units_, function(i) dnorm(out_vec[i], matrix_regression_all[i,], sigma_reg, log=TRUE)) )
    log_weights = log(omega) + log_gauss
    weights = exp(log_weights - apply(log_weights, 1, max))
    
    return(weights)
  } 
  
  
  ####################################################################
  #   ---    parameters    ---
  
  # weights parameters for Intermediates
  Mu_eta_mu = -0.1; Sigma_eta_mu = 0.5
  Gamma1_eta_var = 1; Gamma2_eta_var = 6
  eta_I0_mu = rnorm(dim_W_cov, Mu_eta_mu, Sigma_eta_mu)
  eta_I0_sigma = 1/rgamma(dim_W_cov, shape = Gamma1_eta_var, rate = Gamma2_eta_var) 
  eta_Interm_0 = sapply(1:L_max, function(l) rnorm(dim_W_cov, eta_I0_mu, eta_I0_sigma))
  eta_I1_mu = rnorm(dim_W_cov, Mu_eta_mu, Sigma_eta_mu)
  eta_I1_sigma = 1/rgamma(dim_W_cov, shape = Gamma1_eta_var, rate = Gamma2_eta_var) 
  eta_Interm_1 = sapply(1:L_max, function(l) rnorm(dim_W_cov, eta_I1_mu, eta_I1_sigma))
  # weights parameters for Outcome
  eta_O_mu_0 = rnorm(dim_W_cov, Mu_eta_mu, Sigma_eta_mu)
  eta_O_sigma_0 = 1/rgamma(dim_W_cov, shape = Gamma1_eta_var, rate = Gamma2_eta_var) 
  eta_Outcome_0 = sapply(1:L_max, function(l) rnorm(dim_W_cov, eta_O_mu_0, eta_O_sigma_0))
  eta_O_mu_1 = rnorm(dim_W_cov, Mu_eta_mu, Sigma_eta_mu)
  eta_O_sigma_1 = 1/rgamma(dim_W_cov, shape = Gamma1_eta_var, rate = Gamma2_eta_var) 
  eta_Outcome_1 = sapply(1:L_max, function(l) rnorm(dim_W_cov, eta_O_mu_1, eta_O_sigma_1))
  # regression parameters for Intermediate variable
  Mu_alpha_mu = 0; Sigma_alpha_mu = 5
  Gamma1_alpha_var = 4; Gamma2_alpha_var = 2
  alpha_mu0 = rnorm(dim_I.tr_cov, Mu_alpha_mu, Sigma_alpha_mu)
  alpha_sigma0 = 1/rgamma(dim_I.tr_cov, shape = Gamma1_alpha_var, rate = Gamma2_alpha_var) 
  alpha_0 = t(mvtnorm::rmvnorm(L_max, alpha_mu0, alpha_sigma0*diag(dim_I.tr_cov)))
  alpha_mu1 = rnorm(dim_I.tr_cov, Mu_alpha_mu, Sigma_alpha_mu)
  alpha_sigma1 = 1/rgamma(dim_I.tr_cov, shape = Gamma1_alpha_var, rate = Gamma2_alpha_var) 
  alpha_1 = t(mvtnorm::rmvnorm(L_max, alpha_mu1, alpha_sigma1*diag(dim_I.tr_cov)))
  # variance parameter for Intermediate variable
  gamma1_gamma1 = gamma1_gamma2 = 1;
  gamma2_gamma1 = gamma2_gamma2 = 2;
  Gamma1_sigma_I0 = rgamma(1, shape = gamma1_gamma1, rate = gamma2_gamma1)
  Gamma2_sigma_I0 = rgamma(1, shape = gamma1_gamma2, rate = gamma2_gamma2)
  Sigma_Interm_0 = 1/rgamma(L_max, shape = Gamma1_sigma_I0, rate = Gamma2_sigma_I0) 
  Gamma1_sigma_I1 = rgamma(1, shape = gamma1_gamma1, rate = gamma2_gamma1)
  Gamma2_sigma_I1 = rgamma(1, shape = gamma1_gamma2, rate = gamma2_gamma2)
  Sigma_Interm_1 = 1/rgamma(L_max, shape = Gamma1_sigma_I1, rate = Gamma2_sigma_I1) 
  # regression parameters for Outcome variable
  Mu_beta_mu = 0; Sigma_beta_mu = 5
  Gamma1_beta_var = 4; Gamma2_beta_var = 2
  beta_mu_0 = rnorm(dim_O.tr_cov, Mu_beta_mu, Sigma_beta_mu)
  beta_sigma_0 = 1/rgamma(dim_O.tr_cov, shape = Gamma1_beta_var, rate = Gamma2_beta_var) 
  beta_0 = mvtnorm::rmvnorm(L_max, beta_mu_0, beta_sigma_0*diag(dim_O.tr_cov))
  beta_mu_1 = rnorm(dim_O.tr_cov, Mu_beta_mu, Sigma_beta_mu)
  beta_sigma_1 = 1/rgamma(dim_O.tr_cov, shape = Gamma1_beta_var, rate = Gamma2_beta_var) 
  beta_1 = mvtnorm::rmvnorm(L_max, beta_mu_1, beta_sigma_1*diag(dim_O.tr_cov))
  # variance parameter for Outcome variable
  Gamma1_sigma_O_0 = rgamma(1, shape = gamma1_gamma1, rate = gamma2_gamma1)
  Gamma2_sigma_O_0 = rgamma(1, shape = gamma1_gamma2, rate = gamma2_gamma2)
  Sigma_Outcome_0 = 1/rgamma(L_max, shape = Gamma1_sigma_O_0, rate = Gamma2_sigma_O_0) 
  Gamma1_sigma_O_1 = rgamma(1, shape = gamma1_gamma1, rate = gamma2_gamma1)
  Gamma2_sigma_O_1 = rgamma(1, shape = gamma1_gamma2, rate = gamma2_gamma2)
  Sigma_Outcome_1 = 1/rgamma(L_max, shape = Gamma1_sigma_O_1, rate = Gamma2_sigma_O_1) 
  
  # cluster allocation
  # intermediate:
  weights_Interm_0 = stick_breaking_wights(Cov = weight_Cov_0, eta = eta_Interm_0)
  multinom_Interm_0 <- apply(weights_Interm_0, 1, rmultinom, n=1, size =1)
  Interm_alloc_0 <- c((1:L_max)%*%multinom_Interm_0)
  nL_Interm_0 <- rowSums(multinom_Interm_0)
  
  weights_Interm_1 = stick_breaking_wights(Cov = weight_Cov_1, eta = eta_Interm_1)
  multinom_Interm_1 <- apply(weights_Interm_1, 1, rmultinom, n=1, size =1)
  Interm_alloc_1 <- c((1:L_max)%*%multinom_Interm_1)
  nL_Interm_1 <- rowSums(multinom_Interm_1)
  
  # Outcome:
  weights_Outcome_0 = stick_breaking_wights(Cov = weight_Cov_0, eta = eta_Outcome_0)
  multinom_Outcome_0 <- apply(weights_Outcome_0, 1, rmultinom, n=1, size =1)
  Outcome_alloc_0 <- c((1:L_max)%*%multinom_Outcome_0)
  nL_Outcome_0 <- rowSums(multinom_Outcome_0)
  
  weights_Outcome_1 = stick_breaking_wights(Cov = weight_Cov_1, eta = eta_Outcome_1)
  multinom_Outcome_1 <- apply(weights_Outcome_1, 1, rmultinom, n=1, size =1)
  Outcome_alloc_1 <- c((1:L_max)%*%multinom_Outcome_1)
  nL_Outcome_1 <- rowSums(multinom_Outcome_1)
  
  # dividing units given the allocation
  Cov_interm_0_splitL = lapply(1:L_max, function(l) matrix(interm_Cov_0[Interm_alloc_0==l,], ncol= dim_I.tr_cov))
  Cov_interm_1_splitL = lapply(1:L_max, function(l) matrix(interm_Cov_1[Interm_alloc_1==l,], ncol= dim_I.tr_cov))
  Cov_outcome_0_splitL = lapply(1:L_max, function(l) matrix(outcome_Cov_0[Outcome_alloc_0==l,], ncol = dim_O.tr_cov))
  Cov_outcome_1_splitL = lapply(1:L_max, function(l) matrix(outcome_Cov_1[Outcome_alloc_1==l,], ncol = dim_O.tr_cov))
  
  Interm_0_splitL = lapply(1:L_max, function(l) Interm_0[Interm_alloc_0==l])
  Interm_1_splitL = lapply(1:L_max, function(l) Interm_1[Interm_alloc_1==l])
  Outcome_0_splitL = lapply(1:L_max, function(l) Outcome_0[Outcome_alloc_0==l])
  Outcome_1_splitL = lapply(1:L_max, function(l) Outcome_1[Outcome_alloc_1==l])
  
  ####################################################################
  #   ---    quantities to save    ---
  r_kept = R_tot - R_burnin
  P0_post_allImp <- P1_post_allImp <- matrix(NA, nrow= n_units, ncol= r_kept)
  Y0_P0_post_allImp <- Y1_P0_post_allImp <- matrix(NA, nrow= n_units, ncol= r_kept)
  Y0_P1_post_allImp <- Y1_P1_post_allImp <- matrix(NA, nrow= n_units, ncol= r_kept)
  Decomp_Effects_post <- matrix(NA, nrow= 4, ncol= r_kept)
  row.names(Decomp_Effects_post) <- c("total NDE", "pure NDE", "total NIE", "pure NIE")
  Total_Effects_post <- matrix(NA, nrow= 2, ncol= r_kept)
  row.names(Total_Effects_post) <- c("total_NDE + pure_NIE", "pure_NDE + total_NIE")
  
  ####################################################################
  #   ---    GIBBS START HERE !!    ---
  
  for (r in 1:R_tot){
    
    ### 1- update INTERMEDIATE parameters
    # 1.1: regression paramaters
    var_matrix_inverse_0 = lapply(1:L_max, function(l) solve(crossprod(Cov_interm_0_splitL[[l]]) + (1/alpha_sigma0)*diag(dim_I.tr_cov)))
    mean_vector_0 = lapply(1:L_max, function(l) crossprod(Cov_interm_0_splitL[[l]], Interm_0_splitL[[l]]) + alpha_mu0/alpha_sigma0)
    alpha_0 = sapply(1:L_max, function(l) rmvnorm(1, var_matrix_inverse_0[[l]]%*%mean_vector_0[[l]], var_matrix_inverse_0[[l]]))
    
    var_matrix_inverse_1 = lapply(1:L_max, function(l) solve(crossprod(Cov_interm_1_splitL[[l]]) + (1/alpha_sigma1)*diag(dim_I.tr_cov)))
    mean_vector_1 = lapply(1:L_max, function(l) crossprod(Cov_interm_1_splitL[[l]], Interm_1_splitL[[l]]) + alpha_mu1/alpha_sigma1)
    alpha_1 = sapply(1:L_max, function(l) rmvnorm(1, var_matrix_inverse_1[[l]]%*%mean_vector_1[[l]], var_matrix_inverse_1[[l]]))
    
    # 1.2: variance parameters
    gamma2_partial_sigma_0 = sapply(1:L_max, function(l) sum((Interm_0_splitL[[l]] - Cov_interm_0_splitL[[l]]%*%alpha_0[,l])^2))
    Sigma_Interm_0 = 1/rgamma(L_max, shape = Gamma1_sigma_I0 + nL_Interm_0/2, rate = Gamma2_sigma_I0 + gamma2_partial_sigma_0/2)
    
    gamma2_partial_sigma_1 = sapply(1:L_max, function(l) sum((Interm_1_splitL[[l]] - Cov_interm_1_splitL[[l]]%*%alpha_1[,l])^2))
    Sigma_Interm_1 = 1/rgamma(L_max, shape = Gamma1_sigma_I1 + nL_Interm_1/2, rate = Gamma2_sigma_I1 + gamma2_partial_sigma_1/2)
    
    # 1.3: update 'hyperpar.s' of regression par.s
    sd_alpha_0 = rowSds(alpha_0)^2
    var_hyper_0 = 1/(1/Sigma_alpha_mu + L_max/sd_alpha_0)
    alpha_mu0 = rnorm(dim_I.tr_cov, var_hyper_0*(Mu_alpha_mu/Sigma_alpha_mu + rowSums(alpha_0)/sd_alpha_0), sqrt(var_hyper_0))
    
    alpha_sigma0 = 1/rgamma(dim_I.tr_cov, shape = Gamma1_alpha_var + L_max/2, rate = Gamma2_alpha_var + rowSums((alpha_0 - alpha_mu0)^2)/2) 
    
    sd_alpha_1 = rowSds(alpha_1)^2
    var_hyper_1 = 1/(1/Sigma_alpha_mu + L_max/sd_alpha_1)
    alpha_mu1 = rnorm(dim_I.tr_cov, var_hyper_1*(Mu_alpha_mu/Sigma_alpha_mu + rowSums(alpha_1)/sd_alpha_1), sqrt(var_hyper_1))
    
    alpha_sigma1 = 1/rgamma(dim_I.tr_cov, shape = Gamma1_alpha_var + L_max/2, rate = Gamma2_alpha_var + rowSums((alpha_1 - alpha_mu1)^2)/2) 
    
    # 1.4: update 'hyperpar.s' of variance par.s
    # let's fix Gamma1_sigma_I, make Gamma2_sigma_I varaing
    Gamma2_sigma_I0 = rgamma(1, shape = gamma1_gamma2 + L_max * Gamma1_sigma_I0, rate = gamma2_gamma2 + sum(1/Sigma_Interm_0))
    Gamma2_sigma_I1 = rgamma(1, shape = gamma1_gamma2 + L_max * Gamma1_sigma_I1, rate = gamma2_gamma2 + sum(1/Sigma_Interm_1))
    
    ### 2- cluster allocation for INTERMEDIATE
    # 2.1 posterior weights
    omega_Interm_0 = stick_breaking_wights(Cov = weight_Cov_0, eta = eta_Interm_0)
    weights_Interm_0 = weights_estimation_post(omega = omega_Interm_0,
                                               matrix_reg = interm_Cov_0,
                                               par_reg = alpha_0, 
                                               sigma_reg = Sigma_Interm_0, 
                                               out_vec = Interm_0)
    multinom_Interm_0 = apply(weights_Interm_0, 1, rmultinom, n=1, size =1)
    Interm_alloc_0 = c((1:L_max)%*%multinom_Interm_0)
    nL_Interm_0 = rowSums(multinom_Interm_0)
    
    omega_Interm_1 = stick_breaking_wights(Cov = weight_Cov_1, eta = eta_Interm_1)
    weights_Interm_1 = weights_estimation_post(omega = omega_Interm_1,
                                               matrix_reg = interm_Cov_1,
                                               par_reg = alpha_1, 
                                               sigma_reg = Sigma_Interm_1, 
                                               out_vec = Interm_1)
    multinom_Interm_1 = apply(weights_Interm_1, 1, rmultinom, n=1, size =1)
    Interm_alloc_1 = c((1:L_max)%*%multinom_Interm_1)
    nL_Interm_1 = rowSums(multinom_Interm_1)
    
    # 2.2 save divided matrix for computational efficiency
    Cov_interm_0_splitL = lapply(1:L_max, function(l) matrix(interm_Cov_0[Interm_alloc_0==l,], ncol= dim_I.tr_cov))
    Interm_0_splitL = lapply(1:L_max, function(l) Interm_0[Interm_alloc_0==l])
    
    Cov_interm_1_splitL = lapply(1:L_max, function(l) matrix(interm_Cov_1[Interm_alloc_1==l,], ncol= dim_I.tr_cov))
    Interm_1_splitL = lapply(1:L_max, function(l) Interm_1[Interm_alloc_1==l])
    
    ### 3- augmentation scheme for INTERMEDIATE
    Z_interm_0 = Z_sampling(weights = omega_Interm_0, allocation = Interm_alloc_0)
    Z_interm_1 = Z_sampling(weights = omega_Interm_1, allocation = Interm_alloc_1)
    
    ### 4- update weight-INTERMEDIATE parameters
    # 4.1 update parameters
    Cov_tilde0 = lapply(1:L_max, function(l) matrix(weight_Cov_0[which(Interm_alloc_0>=l),], ncol= dim_W_cov))
    Var_eta0 = lapply(1:L_max, function(l) solve(1/eta_I0_sigma*diag(dim_W_cov) + crossprod(Cov_tilde0[[l]])))
    mean_eta_0 = sapply(1:L_max, function(l) eta_I0_mu/eta_I0_sigma + Z_interm_0[l,which(Interm_alloc_0>=l)]%*%Cov_tilde0[[l]])
    
    eta_Interm_0 = sapply(1:L_max, function(l) mvtnorm::rmvnorm(1, Var_eta0[[l]]%*%mean_eta_0[,l], Var_eta0[[l]]))
    
    Cov_tilde1 = lapply(1:L_max, function(l) matrix(weight_Cov_1[which(Interm_alloc_1>=l),], ncol= dim_W_cov))
    Var_eta1 = lapply(1:L_max, function(l) solve(1/eta_I1_sigma*diag(dim_W_cov) + crossprod(Cov_tilde1[[l]])))
    mean_eta_1 = sapply(1:L_max, function(l) eta_I1_mu/eta_I1_sigma + Z_interm_1[l,which(Interm_alloc_1>=l)]%*%Cov_tilde1[[l]])
    
    eta_Interm_1 = sapply(1:L_max, function(l) mvtnorm::rmvnorm(1, Var_eta1[[l]]%*%mean_eta_1[,l], Var_eta1[[l]]))
    
    
    ### 5- imputation Intermediate NOT needed
    
    ### 6- update OUTCOME parameters
    # 6.1: regression paramaters
    var_matrix_inverse_0 = lapply(1:L_max, function(l) solve(crossprod(Cov_outcome_0_splitL[[l]]) + 
                                                               (1/beta_sigma_0 + .Machine$double.eps)*diag(dim_O.tr_cov)))
    mean_vector_0 = lapply(1:L_max, function(l) crossprod(Cov_outcome_0_splitL[[l]], Outcome_0_splitL[[l]]) + beta_mu_0/beta_sigma_0)
    beta_0 = sapply(1:L_max, function(l) rmvnorm(1, var_matrix_inverse_0[[l]]%*%mean_vector_0[[l]], var_matrix_inverse_0[[l]]))
    
    var_matrix_inverse_1 = lapply(1:L_max, function(l) solve(crossprod(Cov_outcome_1_splitL[[l]]) + 
                                                               (1/beta_sigma_1 + .Machine$double.eps)*diag(dim_O.tr_cov)))
    mean_vector_1 = lapply(1:L_max, function(l) crossprod(Cov_outcome_1_splitL[[l]], Outcome_1_splitL[[l]]) + beta_mu_1/beta_sigma_1)
    beta_1 = sapply(1:L_max, function(l) rmvnorm(1, var_matrix_inverse_1[[l]]%*%mean_vector_1[[l]], var_matrix_inverse_1[[l]]))
    
    # 6.2: variance parameters
    gamma2_partial_sigma_0 = sapply(1:L_max, function(l) sum((Outcome_0_splitL[[l]] - Cov_outcome_0_splitL[[l]]%*%beta_0[,l])^2))
    Sigma_Outcome_0 = 1/rgamma(L_max, shape = Gamma1_sigma_O_0 + nL_Outcome_0/2, rate = Gamma2_sigma_O_0 + gamma2_partial_sigma_0/2)
    
    gamma2_partial_sigma_1 = sapply(1:L_max, function(l) sum((Outcome_1_splitL[[l]] - Cov_outcome_1_splitL[[l]]%*%beta_1[,l])^2))
    Sigma_Outcome_1 = 1/rgamma(L_max, shape = Gamma1_sigma_O_1 + nL_Outcome_1/2, rate = Gamma2_sigma_O_1 + gamma2_partial_sigma_1/2)
    
    # 6.3: update 'hyperpar.s' of regression par.s
    sd_beta_0 = rowSds(beta_0)^2
    var_hyper_0 = 1/(1/Sigma_beta_mu + L_max/sd_beta_0)
    beta_mu_0 = rnorm(dim_O.tr_cov, var_hyper_0*(Mu_beta_mu/Sigma_beta_mu + rowSums(beta_0)/sd_beta_0), sqrt(var_hyper_0))
    
    beta_sigma_0 = 1/rgamma(dim_O.tr_cov, shape = Gamma1_beta_var + L_max/2, rate = Gamma2_beta_var + rowSums((beta_0 -beta_mu_0)^2)/2) 
    beta_sigma_0 <- pmin(pmax(beta_sigma_0, 1e-4), 1e4)   #constrains
    
    sd_beta_1 = rowSds(beta_1)^2
    var_hyper_1 = 1/(1/Sigma_beta_mu + L_max/sd_beta_1)
    beta_mu_1 = rnorm(dim_O.tr_cov, var_hyper_1*(Mu_beta_mu/Sigma_beta_mu + rowSums(beta_1)/sd_beta_1), sqrt(var_hyper_1))
    
    beta_sigma_1 = 1/rgamma(dim_O.tr_cov, shape = Gamma1_beta_var + L_max/2, rate = Gamma2_beta_var + rowSums((beta_1 -beta_mu_1)^2)/2) 
    beta_sigma_1 <- pmin(pmax(beta_sigma_1, 1e-4), 1e4)   #constrains
    
    # 6.4: update 'hyperpar.s' of variance par.s
    # let's fix Gamma1_sigma_O, make Gamma2_sigma_O varing
    Gamma2_sigma_O_0 = rgamma(1, shape = gamma1_gamma2 + L_max * Gamma1_sigma_O_0, rate = gamma2_gamma2 + sum(1/Sigma_Outcome_0))
    Gamma2_sigma_O_1 = rgamma(1, shape = gamma1_gamma2 + L_max * Gamma1_sigma_O_1, rate = gamma2_gamma2 + sum(1/Sigma_Outcome_1))
    
    ### 7- cluster allocation for OUTCOME
    # 7.1 posterior weights
    omega_Outcome_0 = stick_breaking_wights(Cov = weight_Cov_0, eta = eta_Outcome_0)
    weights_Outcome_0 = weights_estimation_post(omega = omega_Outcome_0,
                                                matrix_reg = outcome_Cov_0,
                                                par_reg = beta_0, 
                                                sigma_reg = Sigma_Outcome_0, 
                                                out_vec = Outcome_0)
    multinom_Outcome_0 = apply(weights_Outcome_0, 1, rmultinom, n=1, size =1)
    Outcome_alloc_0 = c((1:L_max)%*%multinom_Outcome_0)
    nL_Outcome_0 = rowSums(multinom_Outcome_0)
    
    omega_Outcome_1 = stick_breaking_wights(Cov = weight_Cov_1, eta = eta_Outcome_1)
    weights_Outcome_1 = weights_estimation_post(omega = omega_Outcome_1,
                                                matrix_reg = outcome_Cov_1,
                                                par_reg = beta_1, 
                                                sigma_reg = Sigma_Outcome_1, 
                                                out_vec = Outcome_1)
    multinom_Outcome_1 = apply(weights_Outcome_1, 1, rmultinom, n=1, size =1)
    Outcome_alloc_1 = c((1:L_max)%*%multinom_Outcome_1)
    nL_Outcome_1 = rowSums(multinom_Outcome_1)
    
    # 7.2 save divided matrix for computational efficiency
    Cov_outcome_0_splitL = lapply(1:L_max, function(l) matrix(outcome_Cov_0[Outcome_alloc_0==l,], ncol = dim_O.tr_cov))
    Outcome_0_splitL = lapply(1:L_max, function(l) Outcome_0[Outcome_alloc_0==l])
    Cov_outcome_1_splitL = lapply(1:L_max, function(l) matrix(outcome_Cov_1[Outcome_alloc_1==l,], ncol = dim_O.tr_cov))
    Outcome_1_splitL = lapply(1:L_max, function(l) Outcome_1[Outcome_alloc_1==l])
    
    ### 8- augmentation scheme for OUTCOME
    Z_outcome_0 = Z_sampling(weights = omega_Outcome_0, allocation = Outcome_alloc_0)
    Z_outcome_1 = Z_sampling(weights = omega_Outcome_1, allocation = Outcome_alloc_1)
    
    ### 9- update weight-OUTCOME parameters
    # 9.1 update parameters
    Cov_tilde0 = lapply(1:L_max, function(l) matrix(weight_Cov_0[which(Outcome_alloc_0>=l),], ncol= dim_W_cov))
    Var_eta0 = lapply(1:L_max, function(l) solve(1/eta_O_sigma_0*diag(dim_W_cov) + crossprod(Cov_tilde0[[l]])))
    mean_eta_0 = sapply(1:L_max, function(l) eta_O_mu_0/eta_O_sigma_0 + Z_outcome_0[l,which(Outcome_alloc_0>=l)]%*%Cov_tilde0[[l]])
    eta_Outcome_0 = sapply(1:L_max, function(l) mvtnorm::rmvnorm(1, Var_eta0[[l]]%*%mean_eta_0[,l], Var_eta0[[l]]))
    
    Cov_tilde1 = lapply(1:L_max, function(l) matrix(weight_Cov_1[which(Outcome_alloc_1>=l),], ncol= dim_W_cov))
    Var_eta1 = lapply(1:L_max, function(l) solve(1/eta_O_sigma_1*diag(dim_W_cov) + crossprod(Cov_tilde1[[l]])))
    mean_eta_1 = sapply(1:L_max, function(l) eta_O_mu_1/eta_O_sigma_1 + Z_outcome_1[l,which(Outcome_alloc_1>=l)]%*%Cov_tilde1[[l]])
    eta_Outcome_1 = sapply(1:L_max, function(l) mvtnorm::rmvnorm(1, Var_eta1[[l]]%*%mean_eta_1[,l], Var_eta1[[l]]))
    
    
    if (r > R_burnin){
      
      # 10 - Principal causal effect 
      # 10.1 SAMPLE imputation step
      # Intermediate 
      log_omega_imputed_0all = log(stick_breaking_wights(Cov = weight_Cov, eta = eta_Interm_0))
      adj_weights_0 = exp(log_omega_imputed_0all - apply(log_omega_imputed_0all, 1, max))
      adj_weights_0 = adj_weights_0/rowSums(adj_weights_0)
      
      log_omega_imputed_1all = log(stick_breaking_wights(Cov = weight_Cov, eta = eta_Interm_1))
      adj_weights_1 = exp(log_omega_imputed_0all - apply(log_omega_imputed_1all, 1, max))
      adj_weights_1 = adj_weights_1/rowSums(adj_weights_1)
      
      reg_tr0 = interm_Cov%*%alpha_0
      reg_tr1 = interm_Cov%*%alpha_1
      
      P0_imp = rowSums(reg_tr0*adj_weights_0)
      P1_imp = rowSums(reg_tr1*adj_weights_1)
      
      # outcome
      Cov_Imp_outcome_Tr0_Interm0[, dim_O.tr_cov] = Cov_Imp_outcome_Tr1_Interm0[, dim_O.tr_cov] = P0_imp
      Cov_Imp_outcome_Tr0_Interm1[, dim_O.tr_cov] = Cov_Imp_outcome_Tr1_Interm1[, dim_O.tr_cov] = P1_imp
      
      log_omega_weights_0all = log(stick_breaking_wights(Cov = weight_Cov, eta = eta_Outcome_0))
      adj_weights_0 = exp(log_omega_weights_0all - apply(log_omega_weights_0all, 1, max))
      adj_weights_0 = adj_weights_0/rowSums(adj_weights_0)
      
      log_omega_weights_1all = log(stick_breaking_wights(Cov = weight_Cov, eta = eta_Outcome_1))
      adj_weights_1 = exp(log_omega_weights_1all - apply(log_omega_weights_1all, 1, max))
      adj_weights_1 = adj_weights_1/rowSums(adj_weights_1)
      
      reg_tr0_i0 = Cov_Imp_outcome_Tr0_Interm0%*%beta_0
      reg_tr0_i1 = Cov_Imp_outcome_Tr0_Interm1%*%beta_0
      reg_tr1_i0 = Cov_Imp_outcome_Tr1_Interm0%*%beta_1
      reg_tr1_i1 = Cov_Imp_outcome_Tr1_Interm1%*%beta_1
      
      Y0_P0_imp = rowSums(reg_tr0_i0*adj_weights_0)
      Y0_P1_imp = rowSums(reg_tr0_i1*adj_weights_0)
      Y1_P0_imp = rowSums(reg_tr1_i0*adj_weights_1)
      Y1_P1_imp = rowSums(reg_tr1_i1*adj_weights_1)
      
      ### Z- save data
      P0_post_allImp[,r-R_burnin] = P0_imp
      P1_post_allImp[,r-R_burnin] = P1_imp
      Y0_P0_post_allImp[,r-R_burnin] = Y0_P0_imp
      Y0_P1_post_allImp[,r-R_burnin] = Y0_P1_imp
      Y1_P0_post_allImp[,r-R_burnin] = Y1_P0_imp
      Y1_P1_post_allImp[,r-R_burnin] = Y1_P1_imp
      
      totNDE = Y1_P1_imp - Y0_P1_imp
      pureNDE = Y1_P0_imp - Y0_P0_imp
      totNIE = Y1_P1_imp - Y1_P0_imp
      pureNIE = Y0_P1_imp - Y0_P0_imp
      
      # 10.2.a SAMPLE Principal causal effect 
      if (method == "sample"){
        Decomp_Effects_post[1,r-R_burnin] = mean(totNDE)
        Decomp_Effects_post[2,r-R_burnin] = mean(pureNDE)
        Decomp_Effects_post[3,r-R_burnin] = mean(totNIE)
        Decomp_Effects_post[4,r-R_burnin] = mean(pureNIE)
        
        Total_Effects_post[1,r-R_burnin] = mean(totNDE + pureNIE)
        Total_Effects_post[2,r-R_burnin] = mean(pureNDE + totNIE)
      }
      
      # 10.2.b POPULATION Principal causal effect 
      if (method == "population") {
        bootstrap_units <- sample(1:n_units, sample_population, replace=TRUE)
        
        Decomp_Effects_post[1,r-R_burnin] = mean(totNDE[bootstrap_units])
        Decomp_Effects_post[2,r-R_burnin] = mean(pureNDE[bootstrap_units])
        Decomp_Effects_post[3,r-R_burnin] = mean(totNIE[bootstrap_units])
        Decomp_Effects_post[4,r-R_burnin] = mean(pureNIE[bootstrap_units])
        
        Total_Effects_post[1,r-R_burnin] = mean(totNDE[bootstrap_units] + pureNIE[bootstrap_units])
        Total_Effects_post[2,r-R_burnin] = mean(pureNDE[bootstrap_units] + totNIE[bootstrap_units])
      }
      
    }
    
    if((r*10) %% R_tot == 0){
      print(paste0(r, "/", R_tot, "[", r/R_tot*100, "%]"))
    }
    
  }
  
  
  ####################################################################
  #   ---    the end    ---
  return(list(Interm_Tr0 = rowMeans(P0_post_allImp), Interm_Tr1 = rowMeans(P1_post_allImp),
              Outcome_Tr0_Interm_Tr0 = rowMeans(Y0_P0_post_allImp), Outcome_Tr0_Interm_Tr1 = rowMeans(Y0_P1_post_allImp),
              Outcome_Tr1_Interm_Tr0 = rowMeans(Y1_P0_post_allImp), Outcome_Tr1_Interm_Tr1 = rowMeans(Y1_P1_post_allImp),
              Direct_Indirect_Effects = Decomp_Effects_post, Total_Effect = Total_Effects_post))
}
