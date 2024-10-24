#library
library(mvtnorm)
library(CholWishart)
library(parallel)
library(truncnorm)
library(invgamma)
library(BNPmix)

SLM_GIbbs<-function(X, Tr, P, Y, R=5000, burnin=1000, n_cluster=10){
  
  # set seed for riproducibility
  set.seed(1)
  
  # ------   preparing variables   ------
  
  # main variables
  T=Tr
  n=length(T)
  
  dD=matrix(rep(0,n*2),ncol=2)
  dD[T==1,2]<-P[T==1]   
  dD[T==0,1]<-P[T==0] 
  
  ########## add X !!!!
  XyC=matrix(c(rep(1,n),dD[,1]),ncol=2)
  XyT=matrix(c(rep(1,n),dD[,1],dD[,2],dD[,1]*dD[,2]),ncol=4)
  
  dD_overall_average = 0*dD
  
  N_XyC = dim(XyC)[2]
  N_XyT = dim(XyT)[2]
  
  beta_yC=matrix(rep(1,N_XyC*R), ncol=R)
  #beta_yC_0=rep(0,N_XyC-1)
  beta_yC_0=rep(0,N_XyC)
  s2_C=10^2
  
  beta_y_d=t(rep(1,R))
  beta_y_d_0=0
  s2_y_d=10^2
  
  beta_yT=matrix(rep(1,N_XyT*R), ncol=R)
  #beta_yT_0=rep(0,N_XyT-1)
  beta_yT_0=rep(0,N_XyT)
  s2_T=10^2
  
  gamma_0 = rep(1,R)
  gamma_0_mu = 3
  gamma_0_s2 = 1
  gamma_1 = rep(1,R)
  gamma_1_mu = 0
  gamma_1_s2 = 1
  s2_3 = 0.5^2
  
  #DP STUFF
  J = n_cluster
  mu_j = matrix(rep(1,J*2*R), nrow=J)             
  G_0_mean = rep(0.5,2)         
  G_0_COV = 0.25^2*diag(2)  
  alpha = 1:R
  
  #SIGMA
  SIGMAdet = rep(0,R)
  SIGMA_m = 3
  SIGMA_0 = diag(2)
  
  pi=matrix(rep(1,J*R),nrow=J)
  
  #INITIAL VALUES:
  beta_yC[,1]=1
  beta_y_d[,1]=1
  beta_yT[,1]=1
  gamma_0[1]=3 
  gamma_1[1]=1 
  SIGMA = matrix(rep(diag(2),J), ncol=2, byrow = TRUE)
  mu_j[,1:2] = matrix(runif(J*2),ncol=2)
  
  pi_star=rbeta(J,1,alpha[1])
  pi_star[J]=1
  pi[1,1]=pi_star[1]
  pi[2:J,1]=sapply(2:J, function(x) pi_star[x]*prod(1-pi_star[1:(x-1)]))
  
  Z=sample(1:J, n, replace=TRUE)
  n_j=1:J
  n_T=sum(T)
  dD[T==0,2]=runif(n-n_T,0,1)
  dD[T==1,1]=runif(n_T,0,1)
  
  component_assign_probs=matrix(rep(Y,J), ncol=J)
  Deviance_Thing=1:R
  
  # -----
  bivnorm<-function(dati,m,sigma){
    sx=sigma[1,1]
    sy=sigma[2,2]
    p=sigma[1,2]/(sx*sy) #error!
    #print(round(c(p,sigma[1,2],sx,sy),4))
    if (abs(p)>=1){
      p=0.99*(p>1)+0.01*(p<0)
    }
    prob_bn=1/(2*3.14159*sx*sy*sqrt(1-p^2))*exp((-1/(2*(1-p^2)))*((dati[1]-m[1])^2/sx^2-2*p*((dati[1]-m[1])/sx)*((dati[2]-m[2])/sy)+(dati[2]-m[2])^2/sy^2))
    return(prob_bn)
  }
  
  # ----
  
  # save information:
  z_late_var=matrix(NA,nrow=n, ncol=R)
  z_late_var[,1]=Z
  post_D_0=matrix(NA,nrow=n, ncol=R)
  post_D_1=matrix(NA,nrow=n, ncol=R)
  post_D_0[,1]=dD[,1]
  post_D_1[,1]=dD[,2]
  post_Y_0=matrix(NA,nrow=n, ncol=R)
  post_Y_1=matrix(NA,nrow=n, ncol=R)
  
  # ---- gibbs sampler: ----
  for (g in 2:R){
    d_cond_D_mean = mu_j[Z[T==1],2*(g-2)+1]+(SIGMA[Z[T==1]*2-1,2]/SIGMA[Z[T==1]*2,2])*
      (dD[T==1,2] - mu_j[Z[T==1],2*(g-2)+2])
    d_cond_D_var = SIGMA[Z[T==1]*2-1,1]-(SIGMA[Z[T==1]*2-1,2]^2/SIGMA[Z[T==1]*2,2])
    
    d_in_y_mean = (Y[T==1]-XyT[T==1, c(1,3)]%*%beta_yT[c(1,3),g-1])/ (matrix(c(rep(1,n_T),XyT[T==1,3]), ncol=2)%*%beta_yT[c(2,4),g-1])
    d_in_y_var = exp(gamma_0[g-1]+XyT[T==1, 3]*gamma_1[g-1])/ (matrix(c(rep(1,n_T),XyT[T==1,3]), ncol=2)%*%beta_yT[c(2,4),g-1])^2 
    
    d_var = 1/(1/d_in_y_var + rep(1,n_T)/d_cond_D_var)
    d_mean = (d_in_y_mean/d_in_y_var + d_cond_D_mean/d_cond_D_var)*d_var
    
    D_cond_d_mean = mu_j[Z[T==0],2*(g-2)+2]+(SIGMA[Z[T==0]*2-1,2]/SIGMA[Z[T==0]*2-1,1])*
      (dD[T==0,1] - mu_j[Z[T==0],2*(g-2)+1])
    D_cond_d_var = SIGMA[Z[T==0]*2,2]-(SIGMA[Z[T==0]*2-1,2]^2/SIGMA[Z[T==0]*2-1,1])
    
    D_var = D_cond_d_var
    D_mean = D_cond_d_mean
    
    dD[T==1,1]=d_mean+sqrt(d_var)*qnorm(sapply(1:n_T, function(i) runif(1,pnorm(0,d_mean[i],sqrt(d_var[i])),pnorm(1,d_mean[i],sqrt(d_var[i])))))
    dD[T==0,2]=D_mean+sqrt(D_var)*qnorm(sapply(1:(n-n_T), function(i) runif(1,pnorm(0,D_mean[i],sqrt(D_var[i])),pnorm(1,D_mean[i],sqrt(D_var[i])))))
    
    # cluster allocation probability
    tmp=matrix(runif(n*2,0,1),nrow=n)
    dD[is.nan(dD)|is.infinite(dD)] = tmp[is.nan(dD)|is.infinite(dD)]
    
    XyC=matrix(c(rep(1,n),dD[,1]),ncol=2)
    XyT=matrix(c(rep(1,n),dD[,1],dD[,2],dD[,1]*dD[,2]),ncol=4)
    
    pi[1,g]=pi_star[1]
    pi[2:J,g]=sapply(2:J, function(x) pi_star[x]*prod(1-pi_star[1:(x-1)]))
    
    for (i in 1:n)
      component_assign_probs[i,]=sapply(1:J, function(j) bivnorm(dD[i,], mu_j[j,(2*(g-1)-1):(2*(g-1))],SIGMA[(2*(j-1)+1):(2*j),])*pi[j,g])
    
    for (i in which(apply(component_assign_probs,1,sum)==0))
      component_assign_probs[i,]=rep(1/J,J)
    
    component_assign_probs=component_assign_probs/apply(component_assign_probs, 1, sum)
    
    component_assign_probs[,J]=1
    for (i in 1:n)
      component_assign_probs[i,(J-1):1]=sapply((J-1):1, function(j) sum(component_assign_probs[i,1:j]))
    
    tmp=matrix(rep(runif(n),J), ncol=J)
    tmp=(tmp>component_assign_probs)
    tmp=1+tmp*matrix(rep(1:J,n),ncol=J,byrow=TRUE)
    Z=apply(tmp, 1, max)
    
    pi_star=sapply(1:J, function(j) rbeta(1,1+sum(Z==j),alpha[g-1]+sum(Z>j)))
    
    #print(paste0(c_n, " / ", round(pi_star,1)))
    
    # updateing Sigma (var matrix in the clusters of D-model)
    for (j in 1:J){
      keep=(Z==j)
      n_j[j]=sum(keep)
      if (n_j[j]==0){
        mu_j[j,(2*(g-1)-1):(2*(g-1))]=rmvnorm(1,G_0_mean,G_0_COV)
      } else {
        if (n_j[j]==1){
          V=solve(n_j[j]*solve(SIGMA[(2*j-1):(2*j),])+solve(G_0_COV))
          mean_o=V%*%(apply(solve(SIGMA[(2*j-1):(2*j),])%*%dD[keep,], 1, sum)+solve(G_0_COV)%*%G_0_mean)
        } else {
          V=solve(n_j[j]*solve(SIGMA[(2*j-1):(2*j),])+solve(G_0_COV))
          mean_o=V%*%(apply(solve(SIGMA[(2*j-1):(2*j),])%*%t(dD[keep,]), 1, sum)+solve(G_0_COV)%*%G_0_mean)}
        mu_j[j,(2*(g-1)+1):(2*(g-1)+2)]=rmvnorm(1,mean_o,V)
      }
    }
    
    for (j in 1:J){
      tmp=dD[Z==j,]-matrix(rep(mu_j[j,(2*(g-1)+1):(2*g)],n_j[j]),ncol=2,byrow=TRUE)
      SIGMA[(j*2-1):(2*j),]=rInvWishart(1,n_j[j]+SIGMA_m,SIGMA_0+t(tmp)%*%tmp)
    }
    
    # updateing beta Y (regression of the mean in Y-model)
    #phi_b2=diag(1/exp(gamma_0[g-1]+dD[T==1,2]*gamma_1[g-1])) 
    #X_tilde_b2=cbind(rep(1,n_T),dD[T==1,2],dD[T==1,2]*dD[T==1,1])
    #Y_tilde_b2=Y[T==1]-dD[T==1,1]*beta_y_d[,g-1]                                                        
    #V_b2=solve(t(X_tilde_b2)%*%phi_b2%*%X_tilde_b2+diag(3)/s2_T)
    #V_b2=(V_b2+t(V_b2))/2
    #M_b2=V_b2%*%(t(X_tilde_b2)%*%phi_b2%*%Y_tilde_b2+beta_yT_0/s2_T)
    #beta_yT[c(1,3:N_XyT),g]=rmvnorm(1,M_b2,V_b2)
    
    #phi_b1=diag(rep(1/exp(gamma_0[g-1]),n-n_T))    
    #X_tilde_b1=rep(1,n-n_T)
    #Y_tilde_b1=Y[T==0]-XyC[T==0,2]*beta_y_d[,g-1]                                                    
    #V_b1=solve(t(X_tilde_b1)%*%phi_b1%*%X_tilde_b1+1/s2_C)
    #V_b1=(V_b1+t(V_b1))/2
    #M_b1=V_b1%*%(t(X_tilde_b1)%*%phi_b1%*%Y_tilde_b1+beta_yC_0/s2_C)
    #beta_yC[1,g]=rmvnorm(1,M_b1,V_b1)
    
    #phi_a=diag(c(rep(1/exp(gamma_0[g-1]),n-n_T),1/exp(gamma_0[g-1]+dD[T==1,2]*gamma_1[g-1])))    
    #X_tilde_a=matrix(c(dD[T==0,1],dD[T==1,1]),nrow=n)
    #Y_tilde_a=c(Y[T==0]-beta_yC[1,g],Y[T==1]-XyT[T==1,c(1,3:N_XyT)]%*%beta_yT[c(1,3:N_XyT),g])                                                         
    #V_a=solve(t(X_tilde_a)%*%phi_a%*%X_tilde_a+1/s2_y_d)
    #V_a=(V_a+t(V_a))/2
    #M_a=V_a%*%(t(X_tilde_a)%*%phi_a%*%Y_tilde_a+beta_y_d_0/s2_y_d)
    #beta_y_d[,g]=rmvnorm(1,M_a,V_a)
    
    #beta_yC[2,g]=beta_y_d[,g]
    #beta_yT[2,g]=beta_y_d[,g]
    
    
    phi_b2=diag(1/exp(gamma_0[g-1]+dD[T==1,2]*gamma_1[g-1])) 
    X_tilde_b2=cbind(rep(1,n_T),dD[T==1,1],dD[T==1,2],dD[T==1,2]*dD[T==1,1])
    Y_tilde_b2=Y[T==1]                                                        
    V_b2=solve(t(X_tilde_b2)%*%phi_b2%*%X_tilde_b2+diag(4)/s2_T)
    V_b2=(V_b2+t(V_b2))/2
    M_b2=V_b2%*%(t(X_tilde_b2)%*%phi_b2%*%Y_tilde_b2+beta_yT_0/s2_T)
    beta_yT[,g]=rmvnorm(1,M_b2,V_b2)
    
    phi_b1=diag(rep(1/exp(gamma_0[g-1]),n-n_T))    
    X_tilde_b1=cbind(rep(1,n-n_T),dD[T==0,1])
    Y_tilde_b1=Y[T==0]                                                   
    V_b1=solve(t(X_tilde_b1)%*%phi_b1%*%X_tilde_b1+diag(2)/s2_C)
    V_b1=(V_b1+t(V_b1))/2
    M_b1=V_b1%*%(t(X_tilde_b1)%*%phi_b1%*%Y_tilde_b1+beta_yC_0/s2_C)
    beta_yC[,g]=rmvnorm(1,M_b1,V_b1)
    
    
    # updateing gamma_0 and gamma_1 (variance of Y-model)
    gamma_0star=rnorm(1,gamma_0_mu, sqrt(gamma_0_s2)) 
    new=log(c(pnorm(Y[T==0],XyC[T==0,]%*%beta_yC[,g],sqrt(exp(gamma_0star))),
              pnorm(Y[T==1],XyT[T==1,]%*%beta_yT[,g],sqrt(exp(gamma_0star+dD[T==1,2]*gamma_1[g-1])))))
    old=log(c(pnorm(Y[T==0],XyC[T==0,]%*%beta_yC[,g],sqrt(exp(gamma_0[g-1]))),
              pnorm(Y[T==1],XyT[T==1,]%*%beta_yT[,g],sqrt(exp(gamma_0[g-1]+dD[T==1,2]*gamma_1[g-1])))))
    old[old==(-Inf)]<- (-100)
    new[new==(-Inf)]<- (-100)
    gamma_0[g]=ifelse(runif(1)<exp(sum(new-old)),gamma_0star,gamma_0[g-1])
    #print(paste0("gamma 0 :",gamma_0[g-1]," / ",gamma_0star," / ",gamma_0[g]))
    
    gamma_1star = rnorm(1,gamma_1_mu, sqrt(gamma_1_s2))
    new = log(pnorm(Y[T==1],XyT[T==1,]%*%beta_yT[,g],sqrt(exp(gamma_0[g]+dD[T==1,2]*gamma_1star)))) 
    old = log(pnorm(Y[T==1],XyT[T==1,]%*%beta_yT[,g],sqrt(exp(gamma_0[g]+dD[T==1,2]*gamma_1[g-1]))))
    old[old==(-Inf)]<- (-100)
    new[new==(-Inf)]<- (-100)
    gamma_1[g]=ifelse(runif(1)<exp(sum(new-old)),gamma_1star,gamma_1[g-1])
    #print(paste0("gamma 1 :",gamma_1[g-1]," / ",sum(new), "/",sum(old)," / ",gamma_1[g]))
    
    alphaSTAR = rgamma(1,1,1)
    tmp = prod(dbeta(pi_star[1:(J-1)],1,alphaSTAR))/prod(dbeta(pi_star[1:(J-1)],1,alpha[g-1]))
    if (is.na(tmp)){tmp=0.5}
    alpha[g]=ifelse(runif(1)<tmp,alphaSTAR,alpha[g-1])
    
    #print(paste0(g, " / ", round(c(prod(dbeta(pi_star[1:(J-1)],1,alphaSTAR)),
    #                               prod(dbeta(pi_star[1:(J-1)],1,alpha[g-1])),
    #                               pi_star),3)))
    
    dD_overall_average=dD_overall_average+dD/(R-1)
    Deviance_Thing[g]= -2*(sum(log(dnorm(Y[T==1], XyT[T==1,]%*%beta_yT[,g],sqrt(exp(gamma_0[g]+dD[T==1,2]*gamma_1[g])))))+
                             sum(log(dnorm(Y[T==0], XyC[T==0,]%*%beta_yC[,g], sqrt(exp(gamma_0[g]))))))
    
    Y_1=rnorm(n, XyT%*%beta_yT[,g],sqrt(exp(gamma_0[g]+dD[T==1,2]*gamma_1[g])))
    Y_0=rnorm(n, XyC[T==0,]%*%beta_yC[,g], sqrt(exp(gamma_0[g])))
    
    # save information
    z_late_var[,g]=Z
    post_D_0[,g]=dD[,1]
    post_D_1[,g]=dD[,2]
    post_Y_0[,g]=Y_0
    post_Y_1[,g]=Y_1
  }
  
  return(list(
    post_P_0_imp=apply(post_D_0[,(burnin+1):R],1,mean, na.rm=TRUE), 
    post_P_1_imp=apply(post_D_1[,(burnin+1):R],1,mean, na.rm=TRUE),
    post_Y_0_imp=apply(post_Y_0[,(burnin+1):R],1,mean, na.rm=TRUE), 
    post_Y_1_imp=apply(post_Y_1[,(burnin+1):R],1,mean, na.rm=TRUE)
  ))
  
}
