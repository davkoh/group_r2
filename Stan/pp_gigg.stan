
data {
  int<lower=1> N; // number of observations
  vector[N] Y; // observations
  int<lower=0> p; // total number of covariates
  int<lower = 0> G; // number of groups
  int pg[p]; // vector of group sizes (of lenth p) // alternative array[p] int pg;
  matrix[N,p] X; //total covariate matrix
  
  // Hyperpriors for the beta' prior
  vector<lower=0>[p] bg; // controls correlation with group shrinkage
  vector<lower=0>[G] ag; // controls group level sparsity
  
  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior

  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;
}

parameters {
  real alpha; // intercept
  real<lower=0> sigma; // noise
  
  // Parameters 
  real<lower=0> tau2; // global shrinkage parameter
  vector<lower=0>[G] gamma2; // group shrinkage factor
  vector<lower=0>[p] lambda2; // local shrinkage parameter
  
  // Group level terms
  vector[p] z_beta;
  }

transformed parameters {
  vector<lower=0>[p] Sigma;
  vector[p] beta;
  
  for (j in 1:p){
    Sigma[j] = tau2 * gamma2[pg[p]] * lambda2[p];
  }
  
  beta  = z_beta .* sqrt(Sigma);
  
}

model {

  // priors
  alpha ~ normal(alpha_mean, alpha_sd);
  z_beta ~ std_normal();
  tau2 ~ student_t(1, 0, sigma);
  sigma ~ normal(0, sigma_sd);
  gamma2 ~ gamma(ag,1);
  lambda2 ~ inv_gamma(bg,1);
  
}

generated quantities{
  vector[N] logliks;
  vector[N] ypred;
  vector[N] mu_pred;
  real postr22;
  for (i in 1:N){
      logliks[i] = normal_lpdf(Y[i]| alpha + X[i,] * beta, sigma);
      ypred[i] = normal_rng( alpha + X[i,] * beta, sigma);
  }
mu_pred = alpha +  X*beta;
postr22 = variance(mu_pred) / (variance(mu_pred) + sigma^2 ); 
}

