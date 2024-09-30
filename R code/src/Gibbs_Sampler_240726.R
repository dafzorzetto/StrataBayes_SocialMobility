StrataBayes_Gibbs <- function(X, X_w, Tr, P, Y, R=5000, burnin=1000, 
                              sigma_beta=NULL, mu_beta=NULL, 
                              gamma_1=1, gamma_2=1,
                              lambda_eta=NULL, mu_eta=NULL,
                              gamma_y_1=1, gamma_y_2=1,
                              mu_eps = NULL, sigma_eps = NULL,
                              dim_cluster=20) {
  
  ### packages
  library(MASS)
  library(truncnorm)
  library(BNPmix)
  library(mvtnorm)
  
  set.seed(1)
  
  n <- nrow(X)
  
  p <- ncol(X)
  p_w <- ncol(X_w)
  
  X_0 <- X[Tr == 0, ,drop = FALSE]
  X_1 <- X[Tr == 1, ,drop = FALSE]
  X_0_w <- X_w[Tr == 0, ,drop = FALSE]
  X_1_w <- X_w[Tr == 1, ,drop = FALSE]
  P_0 <- P[Tr == 0,drop = FALSE]
  P_1 <- P[Tr == 1,drop = FALSE]
  Y_0 <- Y[Tr == 0,drop = FALSE]
  Y_1 <- Y[Tr == 1,drop = FALSE]
  
  n1 <- nrow(X_1); n0 <- nrow(X_0)
  
  
  ### hyperparameters ###
  
  if(is.null(sigma_beta)){
    sigma_beta = diag(2, p, p)
  }
  
  if(is.null(mu_beta)){
    mu_beta = matrix(0, p)
  }
  
  if(is.null(lambda_eta)){
    lambda_eta_0 = diag(2, p+1, p+1)
    lambda_eta_1 = diag(2, p+2, p+2)
  }

  if(is.null(mu_eta)){
    mu_eta_0 = matrix(0, p+1)
    mu_eta_1 = matrix(0, p+2)
  }
  
  if(is.null(sigma_eps)){
    sigma_eps_0 = diag(5, p_w, p_w)
    sigma_eps_1 = diag(5, p_w, p_w)
  }

  if(is.null(mu_eps)){
    mu_eps_0 = matrix(0, p_w)
    mu_eps_1 = matrix(0, p_w)
  } 
  
  # Quantities to store the output
  
  #P0_POST <- P1_POST <- matrix(NA, n, R) #should be just the missing P_t
  P0_POST <- matrix(NA, n, R)
  P1_POST  <- matrix(NA, n, R)
  BETA_0_POST <- BETA_1_POST <- matrix(NA, p, R)
  SIGMA_0_POST <- SIGMA_1_POST <- rep(NA, R)
  Y0_POST <- Y1_POST <- matrix(NA, n, R)
  ETA_0_POST <- array(NA, dim=c(dim_cluster, p+1, R)) 
  ETA_1_POST <- array(NA, dim=c(dim_cluster, p+2, R))
  SIGMA_Y_0_POST <- SIGMA_Y_1_POST <- matrix(NA, dim_cluster, R)
  V0_POST <- V1_POST <- matrix(NA, n, R)
  
  
  # Initialize the sampler
  
  # initialize parameters
  
  beta_0=matrix(0,p)
  beta_1=matrix(0,p)
  sigma_0=5
  sigma_1=10
  eta_0 <- matrix(1, ncol=dim_cluster, nrow=p+1, byrow=TRUE)
  eta_1 <- matrix(1, ncol=dim_cluster, nrow=p+2, byrow=TRUE) #make values smaller
  sigma_y_0 <- rep(1,dim_cluster)
  sigma_y_1 <- rep(1,dim_cluster)
  P_0_update <- ifelse(Tr == 1, median(P[Tr == 0]), P)
  P_1_update <- ifelse(Tr == 0, median(P[Tr == 1]), P)
  Y_0_impute <- ifelse(Tr == 1, median(Y[Tr == 0]), Y)
  Y_1_impute <- ifelse(Tr == 0, median(Y[Tr == 1]), Y)
  V_0_impute <- V_1_impute <- sample(1:dim_cluster, n, replace=TRUE)
  V_0 <- V_0_impute[Tr == 0,drop = FALSE]
  V_1 <- V_1_impute[Tr == 1,drop = FALSE]
  eps_0 <- eps_1 <- rmvnorm(dim_cluster-1, mu_eps_0, sigma_eps_0)
  V_0_all <- V_1_all <- rep(NA,n)
  
  #### PROOF
  #V_0=scenario_2[[1]]$parameters$V[Tr == 0]
  #V_1=scenario_2[[1]]$parameters$V[Tr == 1]
  
  # Precompute the fixed quantities
  
  sigma_beta_inv <- solve(sigma_beta)
  inv_lambda_eta_0 = solve(lambda_eta_0)
  inv_lambda_eta_1 = solve(lambda_eta_1)
  inv_sigma_eps_0 = solve(sigma_eps_0)
  inv_sigma_eps_1 = solve(sigma_eps_1)
  cross_X_0 <- crossprod(X_0)
  cross_X_1 <- crossprod(X_1)
  X_0_P_0 <- t(X_0) %*% P_0
  X_1_P_1 <- t(X_1) %*% P_1
  X_Y_0 <- cbind(X_0,P_0)
  
  #### useful functions
  # probit stick-breaking weights 
  lambda_function=function(X_w,eps){
    U_m=cbind(pnorm(X_w%*%t(eps)),1)
    prod_compl_U_m=cbind(1,t(apply(1-U_m[,-dim_cluster],1,cumprod)))
    return(U_m*prod_compl_U_m)
  }
  # gamma for Z variable
  gamma_function=function(lambda){
    compl_sum_lambda=cbind(1,t(1-apply(lambda[,-dim_cluster],1,cumsum)))
    phi_gamma_argument=lambda/compl_sum_lambda
    phi_gamma_argument[is.infinite(phi_gamma_argument)]<-0.999
    phi_gamma_argument[phi_gamma_argument<=0] <- 0.001
    phi_gamma_argument[phi_gamma_argument>=1] <- 0.999
    return(qnorm(phi_gamma_argument))
  }
  
  # point estimate parition
  estimation_partition<-function(V_post){
    
    # using Variation of Information as loss function
    cluster_post_vi = partition.BNPdens(list(clust=t(V_post)),dist = "VI")$partitions[1,]
    cluster_post_binder = partition.BNPdens(list(clust=t(V_post)),dist = "Binder")$partitions[1,]
    
    return(list(cluster_post_vi=cluster_post_vi,cluster_post_binder=cluster_post_binder))
  }
  
  #progress bar
  pb <- txtProgressBar(style=3)
  
  # Run gibbs sampler
  
  for(r in 1:(R + burnin)){
    
    
    # Print status of chain
    
    setTxtProgressBar(pb, r/(R + burnin))
    
    ### --- post. parameters P(0) ----
    
    ## sampling from the posterior of beta_0
    
    #EPS_beta_new and mu_beta_new are the parameters of the posterior for beta
    EPS_beta_0_new = (solve(sigma_beta_inv + cross_X_0/sigma_0))              
    mu_beta_0_new = EPS_beta_0_new%*%(sigma_beta_inv %*% mu_beta + X_0_P_0/sigma_0)
    
    #sample from the posterior of beta
    beta_0 <- mvrnorm(n = 1, mu = mu_beta_0_new, Sigma = EPS_beta_0_new)
    
    #reg_mu_0 is the mean of the posterior of the P(0)
    #reg_mu_0 = beta_0%*%t(X_0)
    reg_mu_0 = X_0%*%beta_0
    
    #posterior gammas
    gamma_1_post_0 = gamma_1 + n0/2
    gamma_2_post_0 = gamma_2 + sum((P_0 - reg_mu_0)^2)/2
    
    #sample from the posterior sigma
    sigma_0 <- 1/rgamma(1, gamma_1_post_0, gamma_2_post_0)
    
    ### --- post. parameters P(1) ----
    ## sampling from the posterior of beta_1
    
    #EPS_beta_new and mu_beta_new are the parameters of the posterior for beta
    EPS_beta_1_new = (solve(sigma_beta_inv + cross_X_1/sigma_1))              
    mu_beta_1_new = EPS_beta_1_new%*%(sigma_beta_inv %*% mu_beta + X_1_P_1/sigma_1)
    
    #sample from the posterior of beta
    beta_1 <- mvrnorm(n = 1, mu = mu_beta_1_new, Sigma = EPS_beta_1_new)
    
    #reg_mu_1 is the mean of the posterior of the P(1)
    #reg_mu_1 = beta_1%*%t(X_1)
    reg_mu_1 = X_1%*%beta_1
    
    ##sampling from the posterior of gamma
    
    #posterior gammas
    gamma_1_post_1 = gamma_1 + n1/2
    gamma_2_post_1 = gamma_2 + sum((P_1 - reg_mu_1)^2)/2
    
    #sample from the posterior sigma
    sigma_1 <- 1/rgamma(1, gamma_1_post_1, gamma_2_post_1)
    
    ### --- Impute P(0) ----
    
    ###we compute the variance term for P(0) for all components of mixture
    sigma_P0_cl <- sapply(1:dim_cluster, function(v) 
      1/(1/sigma_0 + (eta_1[p+1,v]^2)/sigma_y_0[v]))
    reg_mu_1_mis = X_1%*%beta_0
    
    C = Y_1 - X_1%*%(eta_1[1:p,]) - P_1%*%t(eta_1[p+2,])
    mu_P0_cl <- sigma_P0_cl[V_1]*(reg_mu_1_mis/sigma_0 +eta_1[p+1,V_1]*(diag(C[,V_1]))/sigma_y_0[V_1])
    
    P_0_update[Tr == 1] <- rnorm(n1, mu_P0_cl, sd=sqrt(sigma_P0_cl))
    
    ### --- Impute P(1) ----
    # sample from normale distribution of P(1)
    P_1_update[Tr == 0] <- rnorm(n0, X_0%*%beta_1, sd=sqrt(sigma_1))
    
    
    ### --- Impute V(0) ----
    lambda_0_temp <- lambda_function(X_0_w,eps_0)
    #lambda_0_temp[lambda_0_temp == 0] <- 10^-10
    eta_0_X <- X_Y_0%*%eta_0
    lambda_norm_0 <- matrix(dnorm(rep(Y_0,dim_cluster),c(eta_0_X),rep(sigma_y_0, each=n0),
                                  log=TRUE),
                            ncol=dim_cluster)
    V_0_MN <- sapply(1:n0, function(i)
      tryCatch(
        {
          rmultinom(1,1,exp(lambda_norm_0[i,]+log(lambda_0_temp[i,])))
        },
        error = function(e) {
          # Fill with a 1 in the first slot and zeros in the rest if an error occurs
          c(1, rep(0, length(lambda_norm_0[i, ]) - 1))
        }
      )
      
    )
    V_0 <- c(c(1:dim_cluster)%*%V_0_MN)
    
    ### --- Impute V(1) ----
    X_Y_1 <- cbind(X_1,P_0_update[Tr==1],P_1)
    lambda_1_temp <- lambda_function(X_1_w,eps_1)
    #lambda_1_temp[lambda_1_temp == 0] <- 10^-2
    eta_1_X <- X_Y_1%*%eta_1
    lambda_norm_1 <- matrix(dnorm(rep(Y_1,dim_cluster),c(eta_1_X),rep(sigma_y_1, each=n1),log=TRUE),
                            ncol=dim_cluster)
    V_1_MN <- sapply(1:n1, function(i)
      tryCatch(
        {
          rmultinom(1,1,exp(lambda_norm_1[i,]+log(lambda_1_temp[i,])))
        },
        error = function(e) {
          # Fill with a 1 in the first slot and zeros in the rest if an error occurs
          c(1, rep(0, length(lambda_norm_1[i, ]) - 1))
        }
      )
      
      )
    V_1 <- c(c(1:dim_cluster)%*%V_1_MN)
    
    #### PROOF
    #V_0=scenario_2[[1]]$parameters$V[Tr == 0]
    #V_1=scenario_2[[1]]$parameters$V[Tr == 1]
    
    ### --- Augmentation scheme Z(0) ----
    gamma_0_Z <- gamma_function(lambda_0_temp)
    Z_0 <- sapply(1:n0, function(i) c(rtruncnorm(V_0[i]-1,b=0,mean=gamma_0_Z[i,1:(V_0[i]-1)])*(V_0[i]>1), 
                                      rtruncnorm(1,a=0,mean=gamma_0_Z[i,V_0[i]]), rep(NA,dim_cluster-V_0[i])) )
    
    ### --- Augmentation scheme Z(1) ----
    gamma_1_Z <- gamma_function(lambda_1_temp)
    Z_1 <- sapply(1:n1, function(i) c(rtruncnorm(V_1[i]-1,b=0,mean=gamma_1_Z[i,1:(V_1[i]-1)]),   ##*(V_1[i]>1),
                                      rtruncnorm(1,a=0,mean=gamma_1_Z[i,V_1[i]]), rep(NA,dim_cluster-V_1[i])) )
    
    ### --- Posterior Eps(0) ----
    # product X_tilde and Z_tilde
    prod_X_Z_tilde_0 <- lapply(1:(dim_cluster-1), function(m) (Z_0[m,V_0>=m])%*%(X_0_w[V_0>=m,]))
    Sigma_eps_0 <- lapply(1:(dim_cluster-1), function(m) 
      solve(crossprod(matrix(X_0_w[V_0>=m,], ncol=p_w)) + inv_sigma_eps_0))
    eps_0 <- t(sapply(1:(dim_cluster-1), function(m)
      rmvnorm(1, Sigma_eps_0[[m]]%*%(inv_sigma_eps_0%*%mu_eps_0+t(prod_X_Z_tilde_0[[m]])),Sigma_eps_0[[m]])))
    
    ### --- Posterior Eps(1) ----
    prod_X_Z_tilde_1 <- lapply(1:(dim_cluster-1), function(m) (Z_1[m,V_1>=m])%*%(X_1_w[V_1>=m,]))
    Sigma_eps_1 <- lapply(1:(dim_cluster-1), function(m) 
      solve(crossprod(matrix(X_1_w[V_1>=m,], ncol=p_w)) + inv_sigma_eps_1))
    eps_1 <- t(sapply(1:(dim_cluster-1), function(m)
      rmvnorm(1, Sigma_eps_1[[m]]%*%(inv_sigma_eps_1%*%mu_eps_1+t(prod_X_Z_tilde_1[[m]])),Sigma_eps_1[[m]])))
    
    ### --- post. parameters Y(0) ----
    # eta 0
    crossX_Y0_V <- lapply(1:dim_cluster, function(m)
      crossprod(matrix(X_Y_0[V_0==m, ],ncol=p+1)))
    cross_XY0_V <- lapply(1:dim_cluster, function(m)
      t(matrix(X_Y_0[V_0==m, ],ncol=p+1))%*%Y_0[V_0==m])
    Sigma_eta_V_new <- lapply(1:dim_cluster, function(m)
      solve(crossX_Y0_V[[m]]/sigma_y_0[m] + inv_lambda_eta_0 ))
    
    eta_0 <- sapply(1:dim_cluster, function(m)
      rmvnorm(1, Sigma_eta_V_new[[m]]%*%(inv_lambda_eta_0%*%mu_eta_0 + cross_XY0_V[[m]]/sigma_y_0[m]),Sigma_eta_V_new[[m]]))
    
    # sigma_Y_0
    n0_V <- rowSums(V_0_MN)
    par_gamma <- (Y_0 - X_Y_0%*%eta_0)^2
    sigma_y_0 <- sapply(1:dim_cluster, function(m)
      1/rgamma(1, gamma_y_1+n0_V[m]/2, gamma_y_2 + sum(par_gamma[V_0==m,m])/2))
    
    
    ### --- post. parameters Y(1) ----
    
    crossX_Y1_V <- lapply(1:dim_cluster, function(m)
      crossprod(matrix(X_Y_1[V_1==m, ], ncol=p+2)))
    cross_XY1_V <- lapply(1:dim_cluster, function(m)
      t(matrix(X_Y_1[V_1==m, ],ncol=p+2))%*%Y_1[V_1==m])
    Sigma_eta_V_new <- lapply(1:dim_cluster, function(m)
      solve(crossX_Y1_V[[m]]/sigma_y_1[m] + inv_lambda_eta_1 ))
    
    eta_1 <- sapply(1:dim_cluster, function(m)
      rmvnorm(1, Sigma_eta_V_new[[m]]%*%(inv_lambda_eta_1%*%mu_eta_1 + cross_XY1_V[[m]]/sigma_y_1[m]),Sigma_eta_V_new[[m]]))
    
    n1_V <- rowSums(V_1_MN)
    par_gamma <- (Y_1 - X_Y_1%*%eta_1)^2
    sigma_y_1 <- sapply(1:dim_cluster, function(m)
      1/rgamma(1, gamma_y_1+n1_V[m]/2, gamma_y_2 + sum(par_gamma[V_1==m,m])/2))
    
    ### --- Impute Y(0) ----
    
    # impute Tr=0, using observed data for Tr=1
    lambda_0_mis_temp <- lambda_function(X_1_w,eps_0)
    V_0_mis_MN <- sapply(1:n1, function(i) rmultinom(1,1,lambda_0_mis_temp[i,]))
    V_0_mis <- c(c(1:dim_cluster)%*%V_0_mis_MN)
    V_0_all[Tr==0] <- V_0
    V_0_all[Tr==1] <- V_0_mis
    
    X_Y_0_temp <- cbind(X,P_0_update)
    Y_0_impute <- rnorm(n, rowSums(X_Y_0_temp*t(eta_0[,V_0_all])),sqrt(sigma_y_0[V_0_all]))
    
    ### --- Impute Y(1) ----
    # impute Tr=1, using observed data for Tr=0
    lambda_1_mis_temp <- lambda_function(X_0_w,eps_1)
    V_1_mis_MN <- sapply(1:n0, function(i) rmultinom(1,1,lambda_1_mis_temp[i,]))
    V_1_mis <- c(c(1:dim_cluster)%*%V_1_mis_MN)
    V_1_all[Tr==1] <- V_1
    V_1_all[Tr==0] <- V_1_mis 
    
    X_Y_1_temp <- cbind(X,P_0_update, P_1_update)
    Y_1_impute <- rnorm(n, rowSums(X_Y_1_temp*t(eta_1[,V_1_all])),sqrt(sigma_y_1[V_1_all]))
    
    ### --- Store the sampled values ----
    
    if(r > burnin){
      P0_POST[, r - burnin] <- P_0_update
      P1_POST[, r - burnin] <- P_1_update
      BETA_0_POST[, r - burnin] <- beta_0
      BETA_1_POST[, r - burnin] <- beta_1
      SIGMA_0_POST[r - burnin] <- sigma_0
      SIGMA_1_POST[r - burnin] <- sigma_1
      Y0_POST[, r - burnin] <- Y_0_impute
      Y1_POST[, r - burnin] <- Y_1_impute
      ETA_0_POST[,, r - burnin] <- eta_0
      ETA_1_POST[,, r - burnin] <- eta_1
      SIGMA_Y_0_POST[, r - burnin] <- sigma_y_0
      SIGMA_Y_1_POST[, r - burnin] <- sigma_y_1
      V0_POST[, r - burnin] <- V_0_all
      V1_POST[, r - burnin] <- V_1_all
    }
  }
  
  ### --- point estimation of cluster partition ----
  
  estimation_partition_V0 <- estimation_partition(V0_POST)
  print("estimation_partition_V0 done")
  estimation_partition_V1 <- estimation_partition(V1_POST)
  print("estimation_partition_V1 done")
  
  return(list(P0_POST = P0_POST,
              P1_POST = P1_POST,
              BETA_0_POST=BETA_0_POST,
              BETA_1_POST = BETA_1_POST,
              SIGMA_0_POST = SIGMA_0_POST,
              SIGMA_1_POST = SIGMA_1_POST,
              Y0_POST = Y0_POST,
              Y1_POST = Y1_POST,
              ETA_0_POST = ETA_0_POST,
              ETA_1_POST = ETA_1_POST,
              SIGMA_Y_0_POST = SIGMA_Y_0_POST,
              SIGMA_Y_1_POST = SIGMA_Y_1_POST,
              V0_POST = V0_POST,
              V1_POST = V1_POST,
              estimation_partition_V0 = estimation_partition_V0,
              estimation_partition_V1 = estimation_partition_V1
  ))
  
}
