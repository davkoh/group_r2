

make_dgp <- function(sim_cond){
  name <- sim_cond$dgp
  get(paste0(name, "_dgp"))(sim_cond)
  
}

boss_dgp <- function(sim_cond){
  # Data generating procedure specified in boss et al
  
  #--- Extract settings for data generation
  
  seed <- sim_cond$seed
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  n <- sim_cond$n
  ntest <- sim_cond$ntest
  p <- sim_cond$p
  nu <- sim_cond$nu 
  type <- sim_cond$type
  alpha <- sim_cond$alpha
  #sigma <- sim_cond$sigma #check
  rho_in <- sim_cond$rho_in
  rho_out <- sim_cond$rho_out
  type <- sim_cond$type
  R2 <-  sim_cond$R2
  G <-  sim_cond$Gs # number of groups
  #
  #group_cor_type= rep("AR", groups_n) #type of correlation structure 
  group_ps = rep(sim_cond$p/G, G)  #number of covariates inside each group
  
  # how to generate coefficients inside each group
  group_nus = rep(sim_cond$nu ,G)
  # group_gen_coef_functions= rep(sim_params$gen_coef, groups_n) 
  
  # Mean of X
  mux <- array(0,c(p,1))
  
  # Covariance matrix of X
  block <- diag(p/G)+rho_in -diag(p/G)*rho_in
  covx <- kronecker(diag(G),block)
  covx[covx==0] <- rho_out
  
  # Design matrix X
  
  X <- mvrnorm(n, mux, covx)
  
  # TODO: Find a better solution for this
  if (!is.null(seed)) {
    set.seed(abs(seed-100))
  }
  
  Xtest <- mvrnorm(ntest, mux, covx)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  # Generate beta BOSS dgp
  if (type == "dist"){
    beta <- t(cbind(t(rep(0.5,5)),t(rep(1,5)),t(rep(0,p-p/G))))
  } else {
    beta <- array(0,c(p,1))
    beta[seq(1,p,p/G)] <- c(0.5,1,1.5,2,2)
  }
  
  group_ps <- rep(p/G, G)
  # Generate sigma
  sigma <- as.numeric(sqrt((1-R2)/R2 * t(beta)%*%t(X)%*%X%*%beta / n))
  
  #--- Generate data y 
  # same data for all fits
  
  # coefficients
  alpha <- 0
  #beta <- gen_coef_groups(group_params)
  
  # vector of real parameters
  rtheta <- list(beta= beta, 
                 sigma= sigma, 
                 R2= R2)
  
  # Generate y
  y <-  as.numeric(cbind(rep(1,n), X)%*%c(alpha,beta)+rnorm(n,0, sigma))          
  ytest <-as.numeric(cbind(rep(1,ntest), Xtest)%*%c(alpha,beta)+rnorm(ntest,0, sigma))      
  
  
  dgp_list <- list( seed= seed,
                    y = y, 
                    ytest= ytest, 
                    n= n, 
                    ntest = ntest, 
                    p= p, 
                    X= X, 
                    Xtest= Xtest,
                    rtheta= rtheta,
                    alpha = alpha, 
                    beta = beta, 
                    rho_in = rho_in, 
                    rho_out = rho_out, 
                    covx = covx, 
                    groups_n = G, 
                    group_nus= group_nus,
                    group_ps= group_ps,
                    sigma= sigma)
  
  return(dgp_list)
}

general_dgp <- function(sim_cond){
  
  seed <- sim_cond$seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #------ Generate data
  #--- Extract settings for data generation
  
  # TODO: At the moment all groups are same size with the same behavior
  
  n <- sim_cond$n
  ntest <- sim_cond$ntest
  p <- sim_cond$p
  
  nu <- sim_cond$nu #sparsity per group
  type <- sim_cond$type # type per group
  alpha <- sim_cond$alpha # y intercept
  rho <- sim_cond$rho #rho inside per group 
  
  # TODO: Check rho bw groups
  
  groups_n = sim_cond$Gs # number of groups
  group_rhos= rep(sim_cond$rho, groups_n) #rho in each group
  group_cor_type= rep("AR", groups_n) #type of correlation structure 
  #number of covariates inside each group
  group_ps = rep(sim_cond$p/groups_n, groups_n) #same ps in all groups
  
  # how to generate coefficients inside each group
  group_nus = rep(sim_cond$nu ,groups_n)
  group_gen_coef_functions= rep(sim_params$gen_coef, groups_n) 
  
  # list of parameters that will be used in data generating procedure
  group_params <- list(n=n, 
                       ntest= ntest,
                       p=p,
                       sigma= sigma, 
                       seed= seed,
                       groups_n= groups_n,
                       group_nus= group_nus,
                       group_rhos= group_rhos ,
                       group_cor_type= group_cor_type,
                       group_ps= group_ps, 
                       group_gen_coef_functions= group_gen_coef_functions)
  
  
  # Covariance matrix of X
  covx <- get_sigmaX_grouped(group_params)
  
  # sigma should respond to a prespecified R2
  sigma <- sim_cond$sigma
  
  # Design matrix X
  X <- rmvnorm(n, mean= rep(0,p), sigma=covx)
  
  if (!is.null(seed)) {
    set.seed(abs(seed-100))
  }
  
  Xtest <- mvrnorm(ntest, mean= rep(0,p), sigma=covx)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  

  #--- Generate data y and  b
  # same data for all fits
  
  
  # coefficients
  alpha <- 0
  beta <- gen_coef_groups(group_params)
  
  #vector of real parameters
  rtheta <- list(beta= beta, 
                 #sigma= sim_cond$igma, 
                 R2= sim_cond$R2)
  
  # Generate y
  y <-  as.numeric(cbind(rep(1,n), X)%*%c(alpha,beta)+rnorm(n,0, sigma))          
  ytest <-as.numeric(cbind(rep(1,ntest), Xtest)%*%c(alpha,beta)+rnorm(ntest,0, sigma))      
  
  
  
}