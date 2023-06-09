// GIGG prior
data {
  int<lower=1> N; // number of observations
  vector[N] Y; // observations
  int<lower=0> p; // total number of covariates (includes intercept)
  int<lower = 0> G; // number of groups
  int pg[p-1]; // vector of group sizes (of lenth p) // alternative array[p] int pg;
  matrix[N,p] X; //total covariate matrix (includes column of 1s)
  
  real<lower=0> sigma;  // dispersion parameter
  
  //---- test data
  int<lower=1> Ntest;  // total number of observations
  vector[Ntest] Ytest;  // test set
  matrix[Ntest, p] Xtest;  // population-level design matrix including column of 1s
  
  
  // Hyperpriors for the beta' prior
  vector<lower=0>[p] bg; // controls correlation with group shrinkage
  vector<lower=0>[G] ag; // controls group level sparsity
  
  // sigma prior sd
  //real<lower=0> sigma_sd; // sd of sigma prior
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
  //real alpha; // intercept
  real Intercept;  // temporary intercept for centered predictors
  //real<lower=0> sigma; // noise
  
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
    Sigma[j] = tau2 * gamma2[pg[p]] * lambda2[p] ; // TODO: scale by sdx?
  }
  
  beta  = z_beta .* sqrt(Sigma); 
  
}

model {

  // priors
  //alpha ~ normal(0, 10);
  Intercept ~ normal(0,10);
  z_beta ~ std_normal();
  tau2 ~ student_t(1, 0, sigma);
  //sigma ~ normal(0, sd(Y)); 
  gamma2 ~ gamma(ag,1);
  lambda2 ~ inv_gamma(bg,1);
  
  // likelihood
  Yc ~ normal_id_glm(Xc,Intercept,beta,sigma);
}

generated quantities{
 
  real b_Intercept = Intercept - dot_product(means_X, beta);
  
  vector[N] log_lik; 
  real y_tilde[N];
  vector[N] mu_tilde = rep_vector(0.0, N)+Ymean+Intercept +Xc*beta;
  
  vector[Ntest] log_lik_test; 
  real y_tilde_test[Ntest];
  vector[Ntest] mu_tilde_test = rep_vector(0.0, Ntest)+Ymean+Intercept +Xctest*beta; 
  
  //--- R2
  real<lower=0, upper=1> R2 = variance(mu_tilde) / (variance(mu_tilde) + sigma^2 ); 
  
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

