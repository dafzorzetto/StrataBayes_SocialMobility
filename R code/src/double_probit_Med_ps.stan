
data {
  int<lower=1> N;              // number of observations
  int<lower=1> L_p;            // number of mixture components for P
  int<lower=1> L_y;            // number of mixture components for Y
  int<lower=1> Q;              // number of covariates
  
  matrix[N, Q] X;              // covariates (already standardized)
  vector[N] Tr;                // treatment
  vector[N] ps_hat;            // estimated propensity score
  vector[N] P;                 // univariate post-treatment
  vector[N] Y;                 // univariate final outcome
}


parameters {
  // Centered parameterization
  matrix[L_p-1, Q] beta1_p;
  vector[L_p-1] beta0_p;
  
  matrix[L_y-1, Q] beta1_y;
  vector[L_y-1] beta0_y;
  vector[L_y-1] beta_p_y;
  
  // Component parameters (Gaussian mixture) for P
  vector[L_p] alpha0;                   // intercept - component means
  matrix[L_p, Q] alpha1;                // covariates - component means
  vector[L_p] alpha2;                   // treatment - component means
  vector[L_p] alpha_ps;                 // prop score - component means
  vector<lower=0>[L_p] sigma_p;         // component standard deviations
  
  // Component parameters (Gaussian mixture) for Y
  vector[L_y] gamma0;                   // intercept - component means
  matrix[L_y, Q] gamma1;                // covariates - component means
  vector[L_y] gamma2;                   // treatment - component means
  vector[L_y] gamma3;                   // p_obs - component means
  vector[L_y] gamma_ps;                 // prop score - component means
  vector<lower=0>[L_y] sigma_y;         // component standard deviations
  
  // Hyperparameters for component means (alpha in P distr.)
  real mu0_a;
  vector[Q] mu1_a;
  real mu2_a;
  real mu_ps_a;
  real<lower=0> sigma0_a;
  vector<lower=0>[Q] sigma1_a;
  real<lower=0> sigma2_a;
  real<lower=0> sigma_ps_a;
  
  // Hyperparameters for component means (gamma in Y distr.)
  real mu0_g;
  vector[Q] mu1_g;
  real mu2_g;
  real mu3_g;
  real mu_ps_g;
  real<lower=0> sigma0_g;
  vector<lower=0>[Q] sigma1_g;
  real<lower=0> sigma2_g;
  real<lower=0> sigma3_g;
  real<lower=0> sigma_ps_g;
  
}

transformed parameters {
  
  // Compute observation-specific mixture weights
  matrix[N, L_p] weights_p;
  matrix[N, L_y] weights_y;
  
  for (n in 1:N) {
    vector[L_p-1] eta_p = beta0_p + beta1_p * X[n,]';
    vector[L_p-1] V_p = Phi(eta_p);  // probit transformation
    
    vector[L_y-1] eta_y = beta0_y + beta1_y * X[n,]' + beta_p_y * P[n];
    vector[L_y-1] V_y = Phi(eta_y);  
    
    // Stick-breaking construction for P
    weights_p[n, 1] = V_p[1];
    real remaining_p = 1 - V_p[1];
    for (k_p in 2:(L_p-1)) {
      weights_p[n, k_p] = V_p[k_p] * remaining_p;
      remaining_p *= (1 - V_p[k_p]);
    }
    weights_p[n, L_y] = remaining_p;
    
    // Stick-breaking construction for Y
    weights_y[n, 1] = V_y[1];
    real remaining_y = 1 - V_y[1];
    for (k_y in 2:(L_y-1)) {
      weights_y[n, k_y] = V_y[k_y] * remaining_y;
      remaining_y *= (1 - V_y[k_y]);
    }
    weights_y[n, L_y] = remaining_y;
    
  }
}


model {
  // Priors
  beta0_p ~ normal(-0.5,2);
  to_vector(beta1_p)  ~ normal(0, 2);

  beta0_y ~ normal(-0.5,2);
  to_vector(beta1_y)  ~ normal(0, 2); 
  beta_p_y ~ normal(-0.5,2);
  
  // Component parameters - P
  mu0_a ~ normal(mean(P), 2 * sd(P));
  to_vector(mu1_a) ~ normal(0,5);
  mu2_a ~  normal(0,5);
  mu_ps_a ~  normal(0,5);
  sigma0_a ~ gamma(4, 2);
  to_vector(sigma1_a) ~ gamma(4, 2);
  sigma2_a ~ gamma(4, 2);
  sigma_ps_a ~ gamma(4, 2);
  
  alpha0 ~ normal(mu0_a, sigma0_a);
  alpha2 ~ normal(mu2_a, sigma2_a);
  alpha_ps ~ normal(mu_ps_a, sigma_ps_a);
  sigma_p ~ gamma(2, 1);
  
  for (l in 1:L_p){
    alpha1[l,] ~ normal(mu1_a, sigma1_a);
  }
  
  // Component parameters - Y
  mu0_g ~ normal(mean(Y), 2 * sd(Y));
  to_vector(mu1_g) ~ normal(0,5);
  mu2_g ~  normal(0,5);
  mu3_g ~  normal(0,5);
  mu_ps_g ~  normal(0,5);
  sigma0_g ~ gamma(4, 2);
  to_vector(sigma1_g) ~ gamma(4, 2);
  sigma2_g ~ gamma(4, 2);
  sigma3_g ~ gamma(4, 2);
  sigma_ps_g ~ gamma(4, 2);
  
  gamma0 ~ normal(mu0_g, sigma0_g);
  gamma2 ~ normal(mu2_g, sigma2_g);
  gamma3 ~ normal(mu3_g, sigma3_g);
  gamma_ps ~ normal(mu_ps_g, sigma_ps_g);
  sigma_y ~ gamma(2, 1);
  
  for (l in 1:L_y){
    gamma1[l,] ~ normal(mu1_g, sigma1_g);
  }
  
  // Likelihood: mixture of normals
  for (n in 1:N) {
    vector[L_p] lps_p;
    vector[L_y] lps_y;
    
    vector[L_p] cond_mod_0_p;
    vector[L_p] cond_mod_1_p;
    real P0_sim;
    real P1_sim;
    
    // likelihood for P
    for (k_p in 1:L_p) {
      real alpha_linear = alpha0[k_p] + dot_product(alpha1[k_p,], X[n,]);
      alpha_linear += alpha_ps[k_p] * ps_hat[n];
      
      // likelihood
      lps_p[k_p] = log(weights_p[n, k_p]);
      lps_p[k_p] += normal_lpdf(P[n] | alpha_linear + alpha2[k_p] * Tr[n], sigma_p[k_p]);
      
      // sample P(0) and P(1)
      cond_mod_0_p[k_p] = weights_p[n, k_p] * (alpha_linear);
      cond_mod_1_p[k_p] = weights_p[n, k_p] * (alpha_linear + alpha2[k_p]);
    }
    
    target += log_sum_exp(lps_p);
    
    // likelihood for Y
    for (k_y in 1:L_y) {
      real gamma_linear = gamma0[k_y] + dot_product(gamma1[k_y,], X[n,]) + gamma2[k_y] * Tr[n];
      
      gamma_linear += gamma3[k_y] * P[n] + gamma_ps[k_y] * ps_hat[n];
      
      lps_y[k_y] = log(weights_y[n, k_y]) + normal_lpdf(Y[n] | gamma_linear, sigma_y[k_y]);
    }
    target += log_sum_exp(lps_y);
    
  }
}


generated quantities{
  vector[N] cond_mean_p0;
  vector[N] cond_mean_p1;
  vector[N] cond_mean_y0;
  vector[N] cond_mean_y1;
  
  
  for (n in 1:N) {
    vector[L_p] cond_mod_0_p;
    vector[L_p] cond_mod_1_p;
    vector[L_y] cond_mod_0_y;
    vector[L_y] cond_mod_1_y;
    
    // simulation P
    for (k in 1:L_p) {
      real alpha_linear_common = alpha0[k] + dot_product(alpha1[k,], X[n,]);
      alpha_linear_common += alpha_ps[k] * ps_hat[n];
      
      cond_mod_0_p[k] = weights_p[n, k] * (alpha_linear_common);
      cond_mod_1_p[k] = weights_p[n, k] * (alpha_linear_common + alpha2[k]);
    }
    
    cond_mean_p0[n] = sum(cond_mod_0_p); 
    cond_mean_p1[n] = sum(cond_mod_1_p); 
    
    // simulation Y
    for (j in 1:L_y) {
      real gamma_linear_common = gamma0[j] + dot_product(gamma1[j,], X[n,]) + gamma3[j] * P[n];
      
      gamma_linear_common += gamma_ps[j] * ps_hat[n];
      
      cond_mod_0_y[j] = weights_y[n, j] * (gamma_linear_common);
      cond_mod_1_y[j] = weights_y[n, j] * (gamma_linear_common + gamma2[j]);
    }

    cond_mean_y0[n] = sum(cond_mod_0_y); 
    cond_mean_y1[n] = sum(cond_mod_1_y); 
  }
}
