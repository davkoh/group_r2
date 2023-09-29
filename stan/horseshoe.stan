/* /Horseshoe prior Carvalho, Polson, Scott 2010
beta_i | lambda_i, tau \sim N(0, lambda_i^2 tau^2)
lambda_i \sim HalfCauchy(0,1)
kappa_i= 1/(1+lambda_i^2)

lambda: local shrinkage
tau: global shrinkage
kappa : shrinkage factor for beta

*/
data {
  int<lower=1> N; // Number of observations
  int<lower=1> p; // Number of covariates (includes intercept)
  matrix[N, p] X; // Includes a column of 1s for intercept
  vector[N] y;
  real<lower=0> sigma; // value for sigma
  //---- test data
  int<lower=1> Ntest; // total number of observations
  vector[Ntest] ytest; // test set
  matrix[Ntest, p] Xtest; // population-level design matrix including column of 1s
  
  
  
  int prior_only;  // should the likelihood be ignored?
}

transformed data {
  int pc = p - 1;
  matrix[N, pc] Xc; // centered version of X without an intercept
  vector[pc] means_X; // column means of X before centering
  vector[pc] sds_X; // sds of X before centering
  vector[N] yc;
  real ymean;
  
  matrix[Ntest, pc] Xctest; // centered version of X without an intercept
  vector[pc] means_Xtest; // column means of X before centering
  
  for (i in 2 : p) {
    means_X[i - 1] = mean(X[ : , i]);
    sds_X[i - 1] = sd(X[ : , i]);
    Xc[ : , i - 1] = X[ : , i] - means_X[i - 1];
  }
  
  ymean = mean(y);
  for (i in 1 : N) {
    yc[i] = y[i] - mean(y);
  }
  
  for (i in 2 : p) {
    means_Xtest[i - 1] = mean(Xtest[ : , i]);
    Xctest[ : , i - 1] = Xtest[ : , i] - means_Xtest[i - 1];
  }
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector[pc] zbeta;
  
  vector<lower=0>[pc] lambda;
  real<lower=0> tau;
  //real<lower=0> sigma; // dispersion parameter
  
}

transformed parameters {
  vector[pc] beta;
  beta=zbeta .* lambda * sigma * tau;
}

model {
  
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(y | Xc, Intercept, beta, sigma);
  }
  
  // priors including constants
  target += std_normal_lpdf(zbeta); // prior for zb
  target += cauchy_lpdf(lambda | 0, tau ) //Half Cauchy prior for lambda
    - 1 * cauchy_lcdf(0 | 0, tau);
  
  target += cauchy_lpdf(tau | 0, sigma ) //Half Cauchy prior for tau
    - 1 * cauchy_lcdf(0 | 0, sigma);
  
   //target += student_t_lpdf(sigma | 3, 0, sd(yc))
     //       - 1 * student_t_lccdf(0 | 3, 0, sd(yc)); // scale with sd(y)
            
  target += normal_lpdf(Intercept | 0, 5);  // Intercept
  
}


generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, beta);
  
  vector[N] log_lik;
  array[N] real y_tilde;
  vector[N] mu_tilde = rep_vector(0.0, N) + ymean + Intercept + Xc * beta;
  
  vector[Ntest] log_lik_test;
  array[Ntest] real y_tilde_test;
  vector[Ntest] mu_tilde_test = rep_vector(0.0, Ntest) + ymean + Intercept
                                + Xctest * beta;
  //--- R2
  real<lower=0, upper=1> R2 = variance(mu_tilde)
                              / (variance(mu_tilde) + sigma ^ 2);
  //---y_tilde calc
  
  for (n in 1 : N) {
    log_lik[n] = normal_lpdf(y[n] | mu_tilde[n], sigma);
    y_tilde[n] = normal_rng(mu_tilde[n], sigma); //copy and paste model (executed once per sample) 
  }
  
  //---y_tilde test calc
  for (n in 1 : Ntest) {
    log_lik_test[n] = normal_lpdf(ytest[n] | mu_tilde_test[n], sigma);
    y_tilde_test[n] = normal_rng(mu_tilde_test[n], sigma); //copy and paste model (executed once per sample) 
  }
}

