library(mvtnorm)
library(matrixStats)
library(truncnorm)
library(gtools)         # Dirichlet 

Double_probit_Mediation <- function(Outcome, Interm, Treat, Covariates, PScore = NA,
                             R_tot, R_burnin, L_max, seed=111, 
                             method = "population"){

  set.seed(seed)
  sample_population = 1500
  
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
  if (anyNA(PScore)){
    print("are you sure to do not include NO propensity score?")
  } else if (length(PScore)!=n_units){
    return(print("check dimension of the propensito score")) 
  } else {
    Covariates = cbind(PScore, Covariates)  # propensity score added as covariates
  }
  if(is.na(epsilon)){
    return(print("you have to include epsilon (parameter to define principal causal effect)")) 
  }
  
  # regression matrices
  weight_Cov = cbind(1, Covariates)
  interm_Cov = cbind(weight_Cov, Treat)
  outcome_Cov = cbind(interm_Cov, Interm) 
  dim_W_cov = ncol(weight_Cov)
  dim_I_cov = ncol(interm_Cov)
  dim_O_cov = ncol(outcome_Cov)
  
  which_Tr_0 = which(Treat==0)
  which_Tr_1 = which(Treat==1)
  interm_Cov_Tr1 = interm_Cov[which_Tr_1,]
  interm_Cov_Tr1[,dim_I_cov] <- 0
  n_units_Tr1 = length(which_Tr_1)
  
  interm_Cov_Imp_Tr0 <- interm_Cov_Imp_Tr1 <- interm_Cov
  outcome_Cov_Imp_Tr0 <- outcome_Cov_Imp_Tr1 <- outcome_Cov
  interm_Cov_Imp_Tr0[, dim_I_cov] = outcome_Cov_Imp_Tr0[, dim_O_cov-2] = 0
  interm_Cov_Imp_Tr1[, dim_I_cov] = outcome_Cov_Imp_Tr1[, dim_O_cov-2] = 1
  
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
    
    lim_inf <- matrix(-Inf, nrow=n_units, ncol=L_max)
    lim_sup <- matrix(0, nrow=n_units, ncol=L_max)
    for(i in 1:n_units){
      all_i = allocation[i]
      seq_sup <- if(all_i == L_max) L_max else (all_i+1):L_max
      lim_inf[i, seq_sup] <- NA
      lim_inf[i, all_i] <- 0
      lim_sup[i, seq_sup] <- NA
      lim_sup[i, all_i] <- Inf
    }
    Z_temp = sapply(1:n_units, function(i) rtruncnorm(L_max, a = lim_inf[i,], b = lim_sup[i,], mean = mean_regression[i,], sd = 1))
    return(Z_temp)
  }
  
  
  ####################################################################
  #   ---    parameters    ---
  
  # weights parameters for Intermediates
  Mu_eta_mu = -0.1; Sigma_eta_mu = 0.5
  Gamma1_eta_var = 1; Gamma2_eta_var = 6
  eta_I_mu = rnorm(dim_W_cov, Mu_eta_mu, Sigma_eta_mu)
  eta_I_sigma = 1/rgamma(dim_W_cov, shape = Gamma1_eta_var, rate = Gamma2_eta_var) 
  eta_Interm = sapply(1:L_max, function(l) rnorm(dim_W_cov, eta_I_mu, eta_I_sigma))
  # weights parameters for Outcome
  eta_O_mu = rnorm(dim_W_cov, Mu_eta_mu, Sigma_eta_mu)
  eta_O_sigma = 1/rgamma(dim_W_cov, shape = Gamma1_eta_var, rate = Gamma2_eta_var) 
  eta_Outcome = sapply(1:L_max, function(l) rnorm(dim_W_cov, eta_O_mu, eta_O_sigma))
  # regression parameters for Intermediate variable
  Mu_alpha_mu = 0; Sigma_alpha_mu = 5
  Gamma1_alpha_var = 4; Gamma2_alpha_var = 2
  alpha_mu = rnorm(dim_I_cov, Mu_alpha_mu, Sigma_alpha_mu)
  alpha_sigma = 1/rgamma(dim_I_cov, shape = Gamma1_alpha_var, rate = Gamma2_alpha_var) 
  alpha = t(mvtnorm::rmvnorm(L_max, alpha_mu, alpha_sigma*diag(dim_I_cov)))
  # variance parameter for Intermediate variable
  gamma1_gamma1 = gamma1_gamma2 = 1;
  gamma2_gamma1 = gamma2_gamma2 = 2;
  Gamma1_sigma_I = rgamma(1, shape = gamma1_gamma1, rate = gamma2_gamma1)
  Gamma2_sigma_I = rgamma(1, shape = gamma1_gamma2, rate = gamma2_gamma2)
  Sigma_Interm = 1/rgamma(L_max, shape = Gamma1_sigma_I, rate = Gamma2_sigma_I) 
  # regression parameters for Outcome variable
  Mu_beta_mu = 0; Sigma_beta_mu = 5
  Gamma1_beta_var = 4; Gamma2_beta_var = 2
  beta_mu = rnorm(dim_O_cov, Mu_beta_mu, Sigma_beta_mu)
  beta_sigma = 1/rgamma(dim_O_cov, shape = Gamma1_beta_var, rate = Gamma2_beta_var) 
  beta = mvtnorm::rmvnorm(L_max, beta_mu, beta_sigma*diag(dim_O_cov))
  # variance parameter for Outcome variable
  Gamma1_sigma_O = rgamma(1, shape = gamma1_gamma1, rate = gamma2_gamma1)
  Gamma2_sigma_O = rgamma(1, shape = gamma1_gamma2, rate = gamma2_gamma2)
  Sigma_Outcome = 1/rgamma(L_max, shape = Gamma1_sigma_O, rate = Gamma2_sigma_O) 
  
  # cluster allocation
  # intermediate:
  weights_Interm = stick_breaking_wights(Cov = weight_Cov, eta = eta_Interm)
  multinom_Interm <- apply(weights_Interm, 1, rmultinom, n=1, size =1)
  Interm_alloc <- c((1:L_max)%*%multinom_Interm)
  nL_Interm <- rowSums(multinom_Interm)
  # Outcome:
  weights_Outcome = stick_breaking_wights(Cov = weight_Cov, eta = eta_Outcome)
  multinom_Outcome <- apply(weights_Outcome, 1, rmultinom, n=1, size =1)
  Outcome_alloc <- c((1:L_max)%*%multinom_Outcome)
  nL_Outcome <- rowSums(multinom_Outcome)
  
  # dividing units given the allocation
  Cov_interm_splitL = lapply(1:L_max, function(l) matrix(interm_Cov[Interm_alloc==l,], ncol= dim_I_cov))
  Cov_outcome_splitL = lapply(1:L_max, function(l) matrix(outcome_Cov[Outcome_alloc==l,], ncol = dim_O_cov))
  Interm_splitL = lapply(1:L_max, function(l) Interm[Interm_alloc==l])
  Outcome_splitL = lapply(1:L_max, function(l) Outcome[Outcome_alloc==l])
  
  ####################################################################
  #   ---    quantities to save    ---
  r_kept = R_tot - R_burnin
  P0_post_allImp <- P1_post_allImp <- matrix(NA, nrow= n_units, ncol= r_kept)
  Y0_post_allImp <- Y1_post_allImp <- matrix(NA, nrow= n_units, ncol= r_kept)
  Decomp_Effects_post <- matrix(NA, nrow= 4, ncol= r_kept)
  row.names(Decomp_Effects_post) <- c("total NDE", "pure NDE", "total NIE", "pure NIE")
  Total_Effects_post <- matrix(NA, nrow= 2, ncol= r_kept)
  row.names(Total_Effects_post) <- c("total_NDE + pure_NIE", "pure_NDE + total_NIE")
  ATE_post <- matrix(NA, nrow= 1, ncol= r_kept)
  
  ####################################################################
  #   ---    GIBBS START HERE !!    ---
  
  for (r in 1:R_tot){
    
    ### 1- update INTERMEDIATE parameters
    # 1.1: regression paramaters
    var_matrix_inverse = lapply(1:L_max, function(l) solve(crossprod(Cov_interm_splitL[[l]]) + (1/alpha_sigma)*diag(dim_I_cov)))
    mean_vector = lapply(1:L_max, function(l) crossprod(Cov_interm_splitL[[l]], Interm_splitL[[l]]) + alpha_mu/alpha_sigma)
    alpha = sapply(1:L_max, function(l) rmvnorm(1, var_matrix_inverse[[l]]%*%mean_vector[[l]], var_matrix_inverse[[l]]))
    
    # 1.2: variance parameters
    gamma2_partial_sigma = sapply(1:L_max, function(l) sum((Interm_splitL[[l]] - Cov_interm_splitL[[l]]%*%alpha[,l])^2))
    Sigma_Interm = 1/rgamma(L_max, shape = Gamma1_sigma_I + nL_Interm/2, rate = Gamma2_sigma_I + gamma2_partial_sigma/2)
    
    # 1.3: update 'hyperpar.s' of regression par.s
    sd_alpha = rowSds(alpha)^2
    var_hyper = 1/(1/Sigma_alpha_mu + L_max/sd_alpha)
    alpha_mu = rnorm(dim_I_cov, var_hyper*(Mu_alpha_mu/Sigma_alpha_mu + rowSums(alpha)/sd_alpha), sqrt(var_hyper))
    
    alpha_sigma = 1/rgamma(dim_I_cov, shape = Gamma1_alpha_var + L_max/2, rate = Gamma2_alpha_var + rowSums((alpha - alpha_mu)^2)/2) 
    
    # 1.4: update 'hyperpar.s' of variance par.s
    # let's fix Gamma1_sigma_I, make Gamma2_sigma_I varaing
    Gamma2_sigma_I = rgamma(1, shape = gamma1_gamma2 + L_max * Gamma1_sigma_I, rate = gamma2_gamma2 + sum(1/Sigma_Interm))
    
    ### 2- cluster allocation for INTERMEDIATE
    # 2.1 posterior weights
    omega_Interm = stick_breaking_wights(Cov = weight_Cov, eta = eta_Interm)
    matrix_regression_all = interm_Cov%*%alpha
    log_gauss = do.call(rbind, lapply(1:n_units, function(i) dnorm(Interm[i], matrix_regression_all[i,], Sigma_Interm, log=TRUE)) )
    log_weights = log(omega_Interm) + log_gauss
    weights_Interm = exp(log_weights - apply(log_weights, 1, max))
    multinom_Interm = apply(weights_Interm, 1, rmultinom, n=1, size =1)
    Interm_alloc = c((1:L_max)%*%multinom_Interm)
    nL_Interm = rowSums(multinom_Interm)
    
    # 2.2 save divided matrix for computational efficiency
    Cov_interm_splitL = lapply(1:L_max, function(l) matrix(interm_Cov[Interm_alloc==l,], ncol= dim_I_cov))
    Interm_splitL = lapply(1:L_max, function(l) Interm[Interm_alloc==l])
    
    ### 3- augmentation scheme for INTERMEDIATE
    Z_interm = Z_sampling(weights = omega_Interm, allocation = Interm_alloc)
    
    ### 4- update weight-INTERMEDIATE parameters
    # 4.1 update parameters
    Cov_tilde = lapply(1:L_max, function(l) matrix(weight_Cov[which(Interm_alloc>=l),], ncol= dim_W_cov))
    Var_eta = lapply(1:L_max, function(l) solve(1/eta_I_sigma*diag(dim_W_cov) + crossprod(Cov_tilde[[l]])))
    mean_eta = sapply(1:L_max, function(l) eta_I_mu/eta_I_sigma + Z_interm[l,which(Interm_alloc>=l)]%*%Cov_tilde[[l]])
    
    eta_Interm = sapply(1:L_max, function(l) mvtnorm::rmvnorm(1, Var_eta[[l]]%*%mean_eta[,l], Var_eta[[l]]))
    
    ### 5- imputation Intermediate NOT needed
    
    ### 6- update OUTCOME parameters
    # 6.1: regression paramaters
    var_matrix_inverse = lapply(1:L_max, function(l) solve(crossprod(Cov_outcome_splitL[[l]]) + (1/beta_sigma)*diag(dim_O_cov)))
    mean_vector = lapply(1:L_max, function(l) crossprod(Cov_outcome_splitL[[l]], Outcome_splitL[[l]]) + beta_mu/beta_sigma)
    beta = sapply(1:L_max, function(l) rmvnorm(1, var_matrix_inverse[[l]]%*%mean_vector[[l]], var_matrix_inverse[[l]]))
    
    # 6.2: variance parameters
    gamma2_partial_sigma = sapply(1:L_max, function(l) sum((Outcome_splitL[[l]] - Cov_outcome_splitL[[l]]%*%beta[,l])^2))
    Sigma_Outcome = 1/rgamma(L_max, shape = Gamma1_sigma_O + nL_Outcome/2, rate = Gamma2_sigma_O + gamma2_partial_sigma/2)
    
    # 6.3: update 'hyperpar.s' of regression par.s
    sd_beta = rowSds(beta)^2
    var_hyper = 1/(1/Sigma_beta_mu + L_max/sd_beta)
    beta_mu = rnorm(dim_O_cov, var_hyper*(Mu_beta_mu/Sigma_beta_mu + rowSums(beta)/sd_beta), sqrt(var_hyper))
    
    beta_sigma = 1/rgamma(dim_O_cov, shape = Gamma1_beta_var + L_max/2, rate = Gamma2_beta_var + rowSums((beta -beta_mu)^2)/2) 
    
    # 6.4: update 'hyperpar.s' of variance par.s
    # let's fix Gamma1_sigma_O, make Gamma2_sigma_O varing
    Gamma2_sigma_O = rgamma(1, shape = gamma1_gamma2 + L_max * Gamma1_sigma_O, rate = gamma2_gamma2 + sum(1/Sigma_Outcome))
    
    
    ### 7- cluster allocation for OUTCOME
    # 7.1 posterior weights
    omega_Outcome= stick_breaking_wights(Cov = weight_Cov, eta = eta_Outcome)
    matrix_regression_all = outcome_Cov%*%beta
    log_gauss = do.call(rbind, lapply(1:n_units, function(i) dnorm(Outcome[i], matrix_regression_all[i,], Sigma_Outcome, log=TRUE)) )
    log_weights = log(omega_Outcome) + log_gauss
    weights_Outcome = exp(log_weights - apply(log_weights, 1, max))
    multinom_Outcome = apply(weights_Outcome, 1, rmultinom, n=1, size =1)
    Outcome_alloc = c((1:L_max)%*%multinom_Outcome)
    nL_Outcome = rowSums(multinom_Outcome)
    
    # 7.2 save divided matrix for computational efficiency
    Cov_outcome_splitL = lapply(1:L_max, function(l) matrix(outcome_Cov[Outcome_alloc==l,], ncol = dim_O_cov))
    Outcome_splitL = lapply(1:L_max, function(l) Outcome[Outcome_alloc==l])
    
    ### 8- augmentation scheme for OUTCOME
    Z_outcome = Z_sampling(weights = omega_Outcome, allocation = Outcome_alloc)
    
    ### 9- update weight-OUTCOME parameters
    # 9.1 update parameters
    Cov_tilde = lapply(1:L_max, function(l) matrix(weight_Cov[which(Outcome_alloc>=l),], ncol= dim_W_cov))
    Var_eta = lapply(1:L_max, function(l) solve(1/eta_O_sigma*diag(dim_W_cov) + crossprod(Cov_tilde[[l]])))
    mean_eta = sapply(1:L_max, function(l) eta_O_mu/eta_O_sigma + Z_outcome[l,which(Outcome_alloc>=l)]%*%Cov_tilde[[l]])
    
    eta_Outcome = sapply(1:L_max, function(l) mvtnorm::rmvnorm(1, Var_eta[[l]]%*%mean_eta[,l], Var_eta[[l]]))
    
    
    if (r > R_burnin){
      
      # 10 - Causal effect 
      # 10.1 SAMPLE imputation step
      # Intermediate 
      log_weights = log(omega_Interm)
      adj_weights = exp(log_weights - apply(log_weights, 1, max))
      adj_weights = adj_weights/rowSums(adj_weights)
      reg_tr0 = interm_Cov_Imp_Tr0%*%alpha
      reg_tr1 = interm_Cov_Imp_Tr1%*%alpha
      
      P0_imp = rowSums(reg_tr0*adj_weights)
      P1_imp = rowSums(reg_tr1*adj_weights)
      
      # outcome
      outcome_Cov_Imp_Tr0[which_Tr_1,(dim_O_cov-1)] <- outcome_Cov_Imp_Tr1[which_Tr_1,(dim_O_cov-1)] <- Int_imp
      
      log_weights = log(omega_Outcome)
      adj_weights = exp(log_weights - apply(log_weights, 1, max))
      adj_weights = adj_weights/rowSums(adj_weights)
      reg_tr0 = outcome_Cov_Imp_Tr0%*%beta
      reg_tr1 = outcome_Cov_Imp_Tr1%*%beta
      
      Y0_imp = rowSums(reg_tr0*adj_weights)
      Y1_imp = rowSums(reg_tr1*adj_weights)
      
      ### Z- save data
      P0_post_allImp[,r-R_burnin] = P0_imp
      P1_post_allImp[,r-R_burnin] = P1_imp
      Y0_post_allImp[,r-R_burnin] = Y0_imp
      Y1_post_allImp[,r-R_burnin] = Y1_imp
      
      diff_P = P1_imp - P0_imp
      diff_Y = Y1_imp - Y0_imp
      
      # 10.2.a SAMPLE Principal causal effect 
      if (method == "sample"){
        PCE_post[1,r-R_burnin] = mean(diff_Y[which(diff_P<(-epsilon)) ])
        PCE_post[3,r-R_burnin] = mean(diff_Y[which(diff_P>epsilon) ])
        PCE_post[2,r-R_burnin] = mean(diff_Y[which(diff_P>(-epsilon) & diff_P<epsilon) ])
        
        
        Decomp_Effects_post[1,r-R_burnin] =
        
        Total_Effects_post <- matrix(NA, nrow= 2, ncol= r_kept)
        
        
        #row.names(Decomp_Effects_post) <- c("total NDE", "pure NDE", "total NIE", "pure NIE")
        #row.names(Total_Effects_post) <- c("total_NDE + pure_NIE", "pure_NDE + total_NIE")
        
        
        ATE_post[1,r-R_burnin] = mean(diff_Y)
      }
      
      # 10.2.b POPULATION Principal causal effect 
      if (method == "population") {
        bootstrap_units <- sample(1:n_units, sample_population, replace=TRUE)
        diff_P_boots = diff_P[bootstrap_units]
        diff_Y_boots = diff_Y[bootstrap_units]
        PCE_post[1,r-R_burnin] = mean(diff_Y_boots[which(diff_P_boots<(-epsilon)) ])
        PCE_post[3,r-R_burnin] = mean(diff_Y_boots[which(diff_P_boots>epsilon) ])
        PCE_post[2,r-R_burnin] = mean(diff_Y_boots[which(diff_P_boots>(-epsilon) & diff_P_boots<epsilon) ])
        ATE_post[1,r-R_burnin] = mean(diff_Y_boots)
      }
      
    }
    
    if((r*10) %% R_tot == 0){
      print(paste0(r, "/", R_tot, "[", r/R_tot*100, "%]"))
    }
    
  }
  
  
  ####################################################################
  #   ---    the end    ---
  return(list(Interm_Tr0 = rowMeans(P0_post_allImp), Interm_Tr1 = rowMeans(P1_post_allImp),
              Outcome_Tr0 = rowMeans(Y0_post_allImp), Outcome_Tr1 = rowMeans(Y1_post_allImp),
              Prin_Causal_Effects = PCE_post, Average_Treat_Effect = ATE_post))
}

    
    