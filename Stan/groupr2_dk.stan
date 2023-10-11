data {
  int<lower=1> N;  // total number of observations
  vector[N] y;  // response variable
  int<lower=1> p;  // number of population-level effects, includes intercept
  matrix[N, p] X;  // population-level design matrix, includes a column of 1s
  int<lower=2> G; // number of groups
  //vector[p-1] Ig; // indexes groups
  array[G] int pg; // size of each group
//  real<lower=0> sigma;  // dispersion parameter
  
  // concentration vector of the Dirichlet prior
  // TODO: Change names 
  vector<lower=0>[G] R2D2_alpha_groups; 
  vector<lower=0>[p-1] R2D2_groups_alphas;
  
  
   //---- test data
  int<lower=1> Ntest;  // total number of observations
  vector[Ntest] ytest;  // test set
  matrix[Ntest, p] Xtest;  // population-level design matrix including column of 1s
  
  // data for the R2D2 prior
  real<lower=0> R2D2_mean_R2;  // mean of the R2 prior
  real<lower=0> R2D2_prec_R2;  // precision of the R2 prior
  int prior_only;  // should the likelihood be ignored?
}

transformed data {
  int pc = p - 1;
  matrix[N, pc] Xc;  // centered version of X without an intercept
  vector[pc] means_X;  // column means of X before centering
  vector[pc] sds_X; // sds of X before centering
  vector[N] yc;
  real ymean;
  matrix[Ntest, pc] Xctest;  // centered version of X without an intercept
  vector[pc] means_Xtest;  // column means of X before centering
  
  // TODO: Q should we scale X by group or in total? At the moment the 
  // data that comes is already scaled and centered. 
  
  for (i in 2:p) {
    means_X[i - 1] = mean(X[, i]);
    sds_X[i - 1] =  sd(X[, i]);
    //Xc[, i - 1] = (X[, i] - means_X[i - 1])/sds_X[i - 1];
    Xc[, i - 1] = (X[, i] - means_X[i - 1]);
    
  }
  
  ymean= mean(y);
  for (i in 1:N) {
    yc[i]= y[i]-mean(y);
  }
  
  for (i in 2:p) {
    means_Xtest[i - 1] = mean(Xtest[, i]);
    Xctest[, i - 1] = (Xtest[, i] - means_Xtest[i - 1]) ; 
  }

  
}


  parameters {
  // local parameters for the R2D2 prior
  vector[G] zbeta;
  simplex[G] R2D2_phi_groups; // phi for groups
  real<lower=0 , upper=1> R2D2_R2;  // R2 parameter
  vector<lower=0>[pc] gamma; // Simulate Dirichlets. Long vector of gammas
  real Intercept;  // temporary intercept for centered predictors
  // R2D2 shrinkage parameters
  real<lower=0> sigma;  // dispersion parameter
}



transformed parameters {
  
  real<lower=0> tau2 = R2D2_R2 / (1 - R2D2_R2);
  matrix[N,G] mug;
  vector[G] beta;
  vector[pc] w;
  
  {
    int pos = 1;
  for (g in 1:G){
    w[(pg[g]*(g-1)+1):(g*pg[g])] = segment(gamma, pos, pg[g])/sum(segment(gamma, pos, pg[g]));
    mug[,g] = Xc[,(pg[g]*(g-1)+1):(g*pg[g])] * w[(pg[g]*(g-1)+1):(g*pg[g])];
    pos = pos+pg[g];
}
  }
    
  beta  = zbeta .* sqrt(sigma^2 * tau2 * R2D2_phi_groups);
}

model {
  // priors
   target += std_normal_lpdf(zbeta); //normal distribution zbeta
  target += normal_lpdf(Intercept | 0, 5);  // Intercept  
  // long gammma
  target += gamma_lpdf(gamma | R2D2_groups_alphas, 1); // prior over gamma
  // phi_groups
  target += dirichlet_lpdf(R2D2_phi_groups | R2D2_alpha_groups); // phi_groups ~ dir(alpha)
  
  target += beta_lpdf(R2D2_R2 | R2D2_mean_R2 * R2D2_prec_R2, (1 - R2D2_mean_R2) * R2D2_prec_R2);
  
  sigma ~ student_t(3,0,sd(yc));

  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(yc | mug,Intercept,beta,sigma); //Intercept+Xc*beta
  }
  
}

generated quantities {
  // actual population-level intercept
 // real b_Intercept = Intercept - dot_product(means_X, beta);
  
  vector[N] log_lik; 
  real y_tilde[N];
  vector[N] mu_tilde;
  vector[Ntest] mu_tilde_test;
  matrix[Ntest,G] mugf;
  vector[Ntest] log_lik_test;
  vector[Ntest] y_tilde_test;
  real<lower=0, upper=1> R2;
  
  mu_tilde = rep_vector(0.0, N)+ymean+Intercept +mug*beta;
  
  //--- R2
  R2 = variance(mu_tilde) / (variance(mu_tilde) + sigma^2 ); 
  
  //---y_tilde calc
  for (n in 1:N) {
    log_lik[n] =normal_lpdf( y[n] | mu_tilde[n], sigma); 
    y_tilde[n]= normal_rng(mu_tilde[n], sigma);  //copy and paste model (executed once per sample) 
  }
 
  for (g in 1:G){
    mugf[,g] = Xctest[,(pg[g]*(g-1)+1):(g*pg[g])] * w[(pg[g]*(g-1)+1):(g*pg[g])];
}
  mu_tilde_test = rep_vector(0.0, Ntest)+ymean+Intercept +mugf*beta; 
  
  //---y_tilde test calc
  for (n in 1:Ntest) {
    log_lik_test[n] = normal_lpdf( ytest[n] | mu_tilde_test[n], sigma); 
    y_tilde_test[n]= normal_rng(mu_tilde_test[n], sigma);  //copy and paste model (executed once per sample) 
  }
  
}


