functions {
  /* Efficient computation of the horseshoe prior
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   z: standardized population-level coefficients
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slab regularization parameter
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau ^ 2 * lambda2));
    return z .* lambda_tilde * tau;
  }
}

data {
  int<lower=1> N; // number of time points
  vector[N] Y; // observations
  int<lower = 0> p; // number of covariates assumes intercept
  matrix[N,p] X; //covariate matrix
  
  // data for the horseshoe prior
  real<lower=0> hs_df; // local degrees of freedom
  real<lower=0> hs_df_global; // global degrees of freedom
  real<lower=0> hs_df_slab; // slab degrees of freedom
  real<lower=0> hs_scale_slab; // slab prior scale
  real<lower=0, upper= p> p0; // prior guess on the number of active coefficients
  
  real<lower=0> sigma; // dispersion parameter
   // test data
  int<lower=1> Ntest;  // total number of observations
  vector[Ntest] Ytest;  // test set
  matrix[Ntest, p] Xtest;  // population-level design matrix 
  
  int prior_only;  // should the likelihood be ignored?
  }

transformed data {
  int pc = p - 1;
  matrix[N, pc] Xc;  // centered version of X without an intercept
  vector[pc] means_X;  // column means of X before centering
  vector[pc] sds_X; // sds of X before centering
  vector[N] Yc;
  real Ymean;
  matrix[Ntest, pc] Xctest;  // centered version of X without an intercept
  vector[pc] means_Xtest;  // column means of X before centering
  
  
  for (i in 2:p) {
    means_X[i - 1] = mean(X[, i]);
    sds_X[i - 1] =  sd(X[, i]);
    //Xc[, i - 1] = (X[, i] - means_X[i - 1])/sds_X[i - 1];
    Xc[, i - 1] = (X[, i] - means_X[i - 1]);
    
  }
  
  Ymean= mean(Y);
  for (i in 1:N) {
    Yc[i]= Y[i]-mean(Y);
  }
  
  for (i in 2:p) {
    means_Xtest[i - 1] = mean(Xtest[, i]);
    Xctest[, i - 1] = (Xtest[, i] - means_Xtest[i - 1]) ; 
  }

}

parameters {
  //real<lower=0> sigma; // noise
  vector[pc] zbeta;
  real Intercept;
   
  

  // horseshoe shrinkage parameters
  vector<lower=0>[pc] hs_local;
  real<lower=0> hs_global; // global shrinkage parameter
  real<lower=0> hs_slab; // slab regularization parameter
}

transformed parameters {
  vector[pc] beta = horseshoe(zbeta, hs_local, hs_global, hs_scale_slab ^ 2 * hs_slab);
}

model {
  
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(Yc | Xc, Intercept, beta, sigma); //Intercept+Xc*beta
  }
  

  // priors
  target += std_normal_lpdf(zbeta); //normal distribution zbeta
  target += normal_lpdf(Intercept | 0, 5);  // Intercept 
  hs_global ~ student_t(hs_df_global, 0, p0 / (pc - p0) * sigma / sqrt(N));
  hs_slab ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab);
  hs_local ~ student_t(hs_df, 0, 1);
  //sigma ~ student_t(3,0,sd(Yc));
  
}

generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, beta);
  
  vector[N] log_lik; 
  array[N] real y_tilde;
  vector[N] mu_tilde = rep_vector(0.0, N)+Ymean+Intercept +Xc*beta;
  //--- R2
  real<lower=0, upper=1> R2 = variance(mu_tilde) / (variance(mu_tilde) + sigma^2 );
  
  vector[Ntest] log_lik_test; 
  array[Ntest] real y_tilde_test;
  vector[Ntest] mu_tilde_test = rep_vector(0.0, Ntest)+Ymean+Intercept +Xctest*beta; 
  // lambdas
  vector[pc] lambdas =  hs_scale_slab ^ 2 * hs_slab * square(hs_local) ./ ( hs_scale_slab ^ 2 * hs_slab + hs_global^2 * square(hs_local)); //variances of betas
  

  
  //---y_tilde calc
  for (n in 1:N) {
    log_lik[n] =normal_lpdf( Y[n] | mu_tilde[n], sigma); 
    y_tilde[n]= normal_rng(mu_tilde[n], sigma);  //copy and paste model (executed once per sample) 
  }
  
  //---y_tilde test calc
  for (n in 1:Ntest) {
    log_lik_test[n] =normal_lpdf( Ytest[n] | mu_tilde_test[n], sigma); 
    y_tilde_test[n]= normal_rng(mu_tilde_test[n], sigma);  //copy and paste model (executed once per sample) 
  }
}
