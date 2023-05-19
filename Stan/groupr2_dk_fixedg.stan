data {
  int<lower=1> N; // number of observations
  vector[N] Y; // observations
  int<lower=0> p; // total number of covariates
  int<lower = 0> G; // number of groups
  int<lower=1> pg; // number of group members per group
  matrix[N,p] X; //total covariate matrix
  
  // R2 Stuff
  real<lower=0> mean_R2;  // mean of the R2 prior
  real<lower=0> prec_R2;  // precision of the R2 prior
  vector<lower=0,upper=1>[G] cons; // concentration vector of the Dirichlet prior for group components
  vector<lower=0,upper=1>[p] cons_g; // concentration vector of the Dirichlet prior for subcomponents
  
  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior

  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;
}

parameters {
  real alpha; // intercept
  real<lower=0> sigma; // noise
  
  // R2 Stuff
  simplex[G] psi;
  real<lower=0, upper=1> R2;
  
  // Group level terms
  vector[G] z_beta;
  
  // Group-wise simplexes
  simplex[pg] psig[G];
  }

transformed parameters {
  
  real<lower=0> tau2 = R2 / (1 - R2);
  matrix[N,G] mug;
  vector[p] beta;
  
  
  for (g in 1:G){
    mug[,g] = X[,(pg*(g-1)+1):(g*pg)] * (psig[g,]);
}
  
  beta  = z_beta .* sqrt(sigma^2 * tau2 * psi);
  
}

model {
  int pos1;

  // priors
  alpha ~ normal(alpha_mean, alpha_sd);
  z_beta ~ std_normal();
  psi ~ dirichlet(cons);
  
  for (g in 1:G){
  psig[g,] ~ dirichlet(cons_g[(pg*(g-1)+1):(g*pg)]);
  }
  R2 ~ beta_proportion(mean_R2, prec_R2);
  sigma ~ normal(0, sigma_sd);

  // likelihood
  Y ~ normal_id_glm(mug,alpha,beta,sigma);
}

generated quantities{
  vector[N] logliks;
  vector[N] ypred;
  for (i in 1:N){
      logliks[i] = normal_lpdf(Y[i]| alpha + mug[i,] * beta, sigma);
      ypred[i] = normal_rng( alpha + mug[i,] * beta, sigma);
  }
}

