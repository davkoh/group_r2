# Auxiliary functions 

#------ Plot distributions
plot_dist <- function(dist, bounds, pars, xtype = c("c", "d"), xname="", 
                      prefix = c("d", "p", "q"), parnames = NULL, 
                      package = NULL, ...) {
  xtype <- match.arg(xtype)
  prefix <- match.arg(prefix)
  pos <- -1
  if (!is.null(package)) {
    pos <- asNamespace(package)
  }
  dist_fun <- get(paste0(prefix, dist), pos = pos, mode = "function")
  if (xtype == "c") {
    # continuous
    df <- data.frame(x = seq(bounds[1], bounds[2], 0.001))
  } else if (xtype == "d") {
    # discrete
    df <- data.frame(x = bounds[1]:bounds[2])
  }
  if (!is.null(parnames)) {
    parnames <- paste0(parnames, " = ")
  }
  cnames <- rep(NA, length(pars))
  for (i in seq_along(pars)) {
    tmp <- do.call(dist_fun, c(list(df$x), pars[[i]], list(...)))
    cnames[i] <- paste0("$", parnames, pars[[i]], "$", collapse = ", ")
    df[paste0(parnames, pars[[i]], collapse = ", ")] <- tmp
  }
  
  df <- df %>%
    gather("pars", "dens", -x) %>%
    mutate(pars = factor(pars, unique(pars)))
  
  gg <- ggplot(df, aes(x, dens, color = pars))
  if (xtype == "c") {
    gg <- gg + geom_line(size=1.5)
  } else if (xtype == "d") {
    gg <- gg + 
      geom_linerange(aes(ymin=0, ymax=dens), size = 1) +
      geom_line(size = 0.8, linetype = "dotted", alpha = 0.8)
  }
  
  gg <- gg + 
    #scale_colour_manual(values=my_colors,labels = unname(latex2exp::TeX(cnames)))+
    scale_colour_viridis_d(labels = unname(latex2exp::TeX(cnames)))+
    labs(x = xname, y = "", color = "") + 
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "bottom",
      axis.text.x=element_text(size=20),
      legend.text = element_text(size = 20)
    )
  if (prefix == "p") {
    gg <- gg +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      theme(axis.ticks.y = element_line(), 
            axis.text.y = element_text(),
            axis.line.y = element_line())
  } else if (prefix == "q") {
    gg <- gg +
      scale_y_continuous() +
      theme(axis.ticks.y = element_line(), 
            axis.text.y = element_text(),
            axis.line.y = element_line())
  }
  gg
}


#------- Covariance structures

# Get an ar correlation matrix

get_cor_matrix_ar1 <- function (ar, nobs) {
  # ar: autocorrelation
  # nobs: dimension
  out <- array(0, dim = c(NROW(ar), nobs, nobs))
  fac <- 1/(1 - ar^2)
  pow_ar <- as.list(rep(1, nobs + 1))
  for (i in seq_len(nobs)) {
    pow_ar[[i + 1]] <- ar^i
    out[, i, i] <- fac
    for (j in seq_len(i - 1)) {
      out[, i, j] <- fac * pow_ar[[i - j + 1]]
      out[, j, i] <- out[, i, j]
    }
  }
  out
}

get_sigmaX_ar1 <- function(D, rho){
  # Get
  cor=get_cor_matrix_ar1(ar=rho, nobs=D)
  temp=matrix(0,nrow=D, ncol=D)
  for(i in 1:D){
    temp[i,]=cor[,,i]
  }
  temp
}

cor2cov <- function(R,S){
  # Given a correlation matrix R and a vector of standard deviations S, find the
  # corresponding covariance matrix
  sweep(sweep(R,1,S,"*"), 2, S, "*") #fast
  #diag(S) %*% R %*% diag(S) #slow
}

gen_covmatrix <- function(params){
  # params is a list that contains
  # p : dimension
  # type: type of covariance matrix
  # rho: correlation in case it is an Autoregressive Matrix
  # Generate a covariance matrix 
  
  p= params$p # dimension p of the covariance matrix
  type= params$type 
  
  if(type=="AR"){
    rho= params$rho #correlation 
    # Auto regressive order 1 covariance structure
    covx= get_sigmaX_ar1(p, rho)
    
  }else if(type=="TOE"){
    # Toeplitz covariance structure
    # toe vector specifies the first row of the matrix
    covx=toeplitz(x=params$toe)
    
  }else if( type=="IC"){
    #IC: Intra class correlation
    # Consider that rho should be bw -1/(p-1) and 1
    # For simplicity consider 
    # X_gj= Z_g+ Z_gj , Z_g common norm(0,1), Z_gj group specific norm(0,1)
    
    # TODO: check. I dont remember if this is correct, but it seems so 
    corx= rho*matrix(1,p,p)
    diag(corx)= rep(1, p)
    covx=cor2cov(corx, rep( sqrt(2), p))
    
  }else if(type== "GROUPED"){
    # Create a blocked grouped covariance matrix
    # If this is used then params should have other specific parameters by group
    # For instance the number of groups, size of each group, type of covariance of each group
    # see get_sigmaX_grouped function for more details
    covx = get_sigmaX_grouped(params)  
  }
  
  return(covx)
  
}

get_sigma_sim <- function(params){
  # Given R2 and other params, calculate the corresponding value of sigma
  # See formula in R2D2 or R2D2M2 paper.
  
  R2 <- params$R2
  rho <-  params$rho
  type <-  params$type
  p <- params$p
  nu <- params$nu
  
  #Covariance of X
  covx <- gen_covmatrix(params) 
  sigma_x <- diag(covx)
  
  #Covariance of b
  gencoef_function<- params$gencoef_function
  
  sigma_b <- get(paste0(gencoef_function, "_var"))(list(p=p, nu=nu))
  sigma=sqrt(sum(sigma_x*sigma_b)*(1-R2)/R2)
  
  return(sigma)
  
}

get_sigmaX_grouped <- function(params){
  
  # TODO: improve function!
  # Pay attention to the names of the parameters inside the list params
  
  # groups_n: total number of groups
  # groups_rhos: rho of each group
  # correlation type per group
  # group_ps how many covariates inside each group
  
  ngroups= params$groups_n #number of groups
  rhos = params$group_rhos #correlations of each group
  group_cor_type = params$group_cor_type #correlation type per group
  ps = params$group_ps #number of ps per each group
  
  #list of covariance matrices
  cov_list <- vector(mode = "list", length = ngroups)
  
  for(i in 1:ngroups){
    temp_params <- list(rho= rhos[i],
                        p= ps[i],
                        type= group_cor_type[i])
    
    cov_list[[i]] <- gen_covmatrix(temp_params)
  }
  
  # Create a block matrix of covariances
  as.matrix(Matrix::bdiag(cov_list))
}


#--- Generate coefficients


gen_coef_fixedbs <- function(params){
  # generate a p dimensional vector of coefficients with fixed entries
  # params is a list that contains:
  # p: dimension of the vector
  
  p= params$p
  b= c(rep(2,5), rep(0,p-10) , rep(2,2), rep(0,3)) 
  return(b)
}


gen_coef_sparse_norm <- function(params){
  # generate a p dimensional vector of coefficients with fixed entries
  # Advice: use a normal distribution
  # params is a list that contains:
  # p: dimension of the vector
  # seed: self explanatory
  # 1-nu expected propotion of sparsity
  # TODO: should we also specify the sd of b here?
  
  seed= params$seed
  p= params$p
  nu= params$nu #1-nu is the expected proportion of sparsity
  
  set.seed(seed)
  b= rnorm(p, 0, 3)*rbinom(p, 1 , 1-nu)

  
  return(b)
}


gen_coef_sparse_norm_var <- function(params){
  # Calculate variance of the process
  # p: dimension of the vector
  # seed: self explanatory
  # 1-nu expected propotion of sparsity
  
  p= params$p
  nu= params$nu #sparsity of the vector  p ~ ber(1-nu)
  # assume z_i independent of p_i z_i ~ norm
  # b_i= z_i*p_i 
  # To calculate the variance of b_i and the covariance matrix
  # 1. Montecarlo (but there is an error here!)
  # 2. Law of Total Variance (Might be difficult to calculate)
  
  mu_z= 0
  var_z= 1
  mu_p= 1-nu
  var_p= nu*(1-nu)
  
  #TODO: check formula :)
  #diagonal of the covariance matrix
  var_b= (var_z+mu_z^2)*(var_p+mu_p^2)-mu_z^2*mu_p^2
  var_b= rep(var_b, p) #return diagonal of the covariance matrix
  
  return(var_b)
}


gen_coef_groups <- function(params){
  # generate a p dimensional vector of coefficients with fixed entries
  # Advice: use a normal distribution
  # params is a list that contains:
  # ngroups: number of groups
  # ps: how many covariates inside each group. length(ps) is equal to ngroups
  # nus: vector of expected sparsity per group
  # seed: self explanatory
  # group_gen_coef_functions:  coefficient functions used inside each group
  # TODO: should we also specify the sd of b here?
  
  ngroups= params$groups_n #number of groups
  ps = params$group_ps 
  nus= params$group_nus
  
  seed= params$seed #careful here!
  set.seed(seed)
  seeds= sample(1000000000:.Machine$integer.max,
                size = ngroups)
    
  b=c()
  for(i in 1:ngroups){
    p= ps[i] 
    nu= nus[i]
    params_temp= list(p=p, nu=nu, seed= seeds[i])
    #get coefficient function
    gen_coef_fn <- get( params$group_gen_coef_functions[i], 
                        mode = "function")
    b=c(b,gen_coef_fn(params_temp))
  }
  b
}

gen_coef_groups_var <- function(params){
  # Calculate variance of the process
  # ps: how many covariates inside each group. length(ps) is equal to ngroups
  # nus: vector of expected sparsity per group
  # seed: self explanatory
  # group_gen_coef_functions:  coefficient functions used inside each group
  
  #sparsity of the vector  p ~ ber(1-nu)
  # assume z_i independent of p_i z_i ~ norm
  # b_i= z_i*p_i 
  # To calculate the variance of b_i and the covariance matrix
  # 1. Montecarlo (but there is an error here)
  # 2. Law of Total Variance
  
  
  # In any case, 
  
  ngroups= params$groups_n #number of groups
  ps = params$group_ps 
  nus= params$group_nus
  
  var_b=c()
  for(i in 1:ngroups){
    p= ps[i] 
    nu= nus[i]
    params_temp= list(p=p, nu=nu, seed= params$seed)
    
    #calculate variance per group
    gen_coef_fn_var <- get( paste0(params$group_gen_coef_functions[i], "_var"),
                            mode = "function")
    #return vector of variances
    var_b=c(var_b,gen_coef_fn_var(params_temp))
    
  }
  var_b
}

#---------- Simulation Functions


dataset_cond_sim <- function(sim_cond,
                             sim_params,
                             smqoi, 
                             seed= NULL,
                             path= NULL,
                             ncores= 1,
                             special_name=NULL){
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  seed_list <- sample(100000000:.Machine$integer.max,
                      size = sim_params$nsims)
  
  ParallelLogger::clearLoggers()
  ParallelLogger::addDefaultConsoleLogger(name=paste0("dataset_cond_sim_", as.character(seed)))
  ParallelLogger::addDefaultErrorReportLogger(name=paste0("dataset_cond_sim_", as.character(seed)))
  on.exit(unregisterLogger(paste0("dataset_cond_sim_", as.character(seed))))
  
  if (file.exists(paste0(paste(path, paste0(special_name,"_",sim_cond$id), sep = "/"), ".RDS"))) {
    return(readRDS(paste0(paste(path, paste0(special_name,"_",sim_cond$id), sep = "/"), ".RDS")))
  } else {
    if (ncores > 1) {
      # Multiprocessing setup
      
      cluster <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cluster)
      
      parallel::clusterExport(cl= cluster, c('sim_params', 'smqoi'))
      
      parallel::clusterEvalQ(cl = cluster, {
        library(brms)
        library(tidybayes)
        library(cmdstanr)
        library(tibble)
        library(dplyr)
        library(mvtnorm)
        library(LaplacesDemon)
        library(parallel)
        library(Matrix)
        library(gtools)
        library(ParallelLogger)
        
        library(MASS)
        
        #Run important scripts
        source("R/aux_functions.R")
        source("R/dgp_functions.R")
        source("R/big_sim_full_simulation.R")
        source("R/R2D2_alpha_gen.R")
        source("R/big_sim_mcmc_params.R")
        source("R/stan_fits.R")
        
        cmdstanr::set_cmdstan_path(sim_params$cmdstan_path)
        
      })
      
      `%dopar%` <- foreach::`%dopar%`
      
      # Multiprocessing run
      results <- foreach::foreach(
        par_seed = seed_list
      ) %dopar% {
        cond_sim(
          sim_params = sim_params,
          sim_cond = sim_cond,
          smqoi = smqoi,
          seed = par_seed
        )

        
      }
      
      # Multiprocessing teardown
      parallel::stopCluster(cluster)
    } 
    # TODO : add single core processing (someday)
    final_result <- do.call(rbind, results)
    #final_result$data_config_seed <- seed
    
    if (!is.null(path)) {
      saveRDS(final_result,
              paste0(paste(path, paste0(special_name,"_",sim_cond$id), 
                           sep = "/"), ".RDS"))
    }
    return(final_result) 
  }
}

#Simulate data and summarise a fit given a simulation condition
cond_sim <-  function(sim_params, sim_cond, smqoi, seed = NULL){
  
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  #------ Generate data
  
  # Data generating procedure is needed
  # Data gen procedure should be inside sim_cond
  # The function should return a list
  
  sim_cond$seed <- seed #need the seed to control the dgp
  dgp_list <- make_dgp(sim_cond)
  
  #------ Fit different models
  
  fits_params <- sim_params$fits_params 
  nfits <- fits_params$nfits #number of fits
  fits_list <- fits_params$fits_list #fits to considered
  names_list <- fits_params$names_list
  nnames <- fits_params$nnames
  
  
  summary_list <- vector(mode = "list", length = nnames)
  params_list <- vector(mode = "list", length = nnames)
  
  #--- Create a temporary folder to store stan csv files
  
  temp_directory <- paste0("temp_stan_files/",seed)
  dir.create(temp_directory)
  dgp_list$temp_directory <- temp_directory
  
  
  #Cycle through models
  for(i in 1:nnames){
    
    
    mcmc_params <- make_mcmc_params(names_list[i], 
                                    fits_params$mcmc_params[[i]], 
                                    data_gen_params = dgp_list)
    
    mcmc_params$iter <-  fits_params$mcmc_params[[i]]$iter
    
    params_list[[i]] <- append(mcmc_params,
                               dgp_list)
    
    
    fitfn <- get(paste0(fits_list[i],"fit")) # Get fit function
    fit <- fitfn(params_list[[i]]) # Run stan fit
    
    #Summarize
    if(fit$fit.error){
      summary_list[[i]] <- NULL
    }else{
      fit_summary_params <- list(fit= fit$fit,
                                 rtheta= dgp_list$rtheta,
                                 standat= params_list[[i]],
                                 voi=smqoi$voi,
                                 moi= smqoi$moi,
                                 seed= seed,
                                 probsoi= smqoi$probsoi)
      
      summary_list[[i]] <- myfitsummary(fit_summary_params)
    }      
    
  }
  
  #delete the temporary stan files
  unlink(temp_directory, recursive = TRUE)
  names(summary_list) <- names_list
  
  return(summary_list)
  
}

#--- Summary functions

#--Percentile estimator
invq <- function(theta,draws){
  #estimator of rqs, such that given theta
  # P( X <= theta)= rqs
  rqs= c()
  for(j in 1:length(theta)){
    rqs[j]=mean(draws[,,j]<theta[j])
  }
  return(rqs) 
}


#--Calculate posterior rmse observations
prmse <-function(standat, fit){
  #posterior rmse
  y=standat$y
  ytest=standat$ytest
  N=standat$n
  Ntest=standat$ntest
  
  # tildes
  
  tildes=fit$draws(c("y_tilde",
                     "mu_tilde",
                     "y_tilde_test" ,
                     "mu_tilde_test"),
                   format="draws_matrix")
  
  niter = dim(tildes)[1]
  y_tilde = tildes[,1:N]
  mu_tilde = tildes[,(N+1):(2*N)]
  y_tilde_test = tildes[,(2*N+1):(2*N+Ntest)]
  mu_tilde_test = tildes[,(2*N+1+Ntest):(2*N+2*Ntest)]
  
  #train mse
  ytemp=matrix(y, byrow = TRUE, ncol = N, nrow = niter) 
  mse1_train=mean((y-colMeans(y_tilde))^2)
  mse2_train=mean((ytemp-y_tilde)^2)
  mse3_train=mean((y-colMeans(mu_tilde))^2)
  mse4_train=mean((ytemp-mu_tilde)^2)
  
  #test mse
  ytesttemp=matrix(ytest, byrow = TRUE, ncol = Ntest, nrow = niter) 
  mse1_test=mean((ytest-colMeans(y_tilde_test))^2)
  mse2_test=mean((ytesttemp-y_tilde_test)^2)
  mse3_test=mean((ytest-colMeans(mu_tilde_test))^2)
  mse4_test=mean((ytesttemp-mu_tilde_test)^2)
  
  mses_train= c(mse1_train,mse2_train,mse3_train,mse4_train)
  mses_test= c(mse1_test,mse2_test, mse3_test, mse4_test)
  
  return(sqrt(c(mses_train, mses_test)))
  
}  

prmse_theta <- function(fit, theta_name, real_theta) {
  # Calculate posterior RMSE for theta
  theta_hat = as_draws_matrix(fit$draws(variables = theta_name))
  theta_real =matrix( real_theta, 
                      byrow = TRUE, 
                      nrow =dim(theta_hat)[1],
                      ncol = dim(theta_hat)[2])
  
  temp <- rowSums(((theta_hat - theta_real))*((theta_hat - theta_real)))
  sqrt(mean(temp))
  
}

prmse_theta0 <- function(fit, theta_name, real_theta) {
  # Posterior RMSE for null theta
  
  Index <- real_theta==0
  
  if( sum(Index) > 0 ){
    real_theta_0 <- real_theta[Index]
    
    theta_hat_0 <- as_draws_matrix(fit$draws(variables = theta_name))[,Index]
    
    theta_temp= matrix( real_theta_0, 
                        byrow = TRUE, 
                        nrow =dim(theta_hat_0)[1], 
                        ncol = dim(theta_hat_0)[2])
    
    temp <- rowSums(((theta_hat_0 - theta_temp))*((theta_hat_0 - theta_temp)))
    return(sqrt(mean(temp)))
  }else{
    return(NA)
  }
}


prmse_thetan0 <- function(fit, theta_name, real_theta) {
  #posterior RMSE for non null theta
  Index <- real_theta!=0
  if( sum(Index) > 0 ){
    real_theta_n0 <- real_theta[Index]
    
    theta_hat_n0 <- as_draws_matrix(fit$draws(variables = theta_name))[,Index]
    
    theta_temp= matrix( real_theta_n0, 
                        byrow = TRUE, 
                        nrow =dim(theta_hat_n0)[1], 
                        ncol = dim(theta_hat_n0)[2])
    
    temp <- rowSums(((theta_hat_n0 - theta_temp))*((theta_hat_n0 - theta_temp)))
    return(sqrt(mean(temp)))
  }else{
    return(NA)
  }
  
}


prmse_theta_pp <- function(fit, theta_name, real_theta){
  # Per parameter RMSE
  # prmse for each coefficient
  
  # theta_hat_i is the ith posterior draw
  # This is our favorite one c:
  # We want to have comparability across different ps
  temp= c()
  theta_hat = as_draws_matrix(fit$draws(variables = theta_name))
  for(i in 1:length(real_theta)){
    #rmse per parameter
    temp[i]= sqrt(mean((theta_hat[, i]-real_theta[i])^2)) # how far theta_hat_i is from real_theta
  }
  
  # Average of rmses across all posterior samples
  mean(temp)
  
}


prmse_theta_pp0 <- function(fit, theta_name, real_theta){
  # per parameter RMSE on null coefficients only
  # prmse for each coefficient
  
  # theta_hat_i is the ith posterior draw
  # This is our favorite one c:
  # We want to have comparability across different ps
  Index <- real_theta==0
  
  if(sum(Index)>0){
    temp= c()
    theta_hat = as_draws_matrix(fit$draws(variables = theta_name))
    theta_hat = theta_hat[,Index]
    for(i in 1:length(real_theta[Index])){
      #rmse per parameter
      temp[i]= sqrt(mean((theta_hat[, i]-real_theta[i])^2)) # how far theta_hat_i is from real_theta
    }
    # Average of rmses across all posterior samples
    return(mean(temp))
    
  }else{
    return(NA)
  }
  
}

prmse_theta_ppn0 <- function(fit, theta_name, real_theta){
  
  # per parameter RMSE
  # prmse for each coefficient
  
  # This is our favorite one c:
  # We want to have comparability across different ps
  Index <- real_theta!=0
  
  if(sum(Index)>0){
    temp= c()
    theta_hat = as_draws_matrix(fit$draws(variables = theta_name))
    theta_hat = theta_hat[,Index]
    for(i in 1:length(real_theta[Index])){
      #rmse per parameter
      temp[i]= sqrt(mean((theta_hat[, i]-real_theta[i])^2)) # how far theta_hat_i is from real_theta
    }
    # Average of rmses across all posterior samples
    return(mean(temp))
    
  }else{
    return(NA)
  }
}



myfitsummary <- function(fit_summary_params){
  #fit: stan fit
  #qoi: quantities of interest
  #variables of interest:
  
  fit = fit_summary_params$fit
  rtheta = fit_summary_params$rtheta
  standat =fit_summary_params$standat
  p = standat$p
  voi =fit_summary_params$voi
  moi = fit_summary_params$moi 
  probsoi = fit_summary_params$probsoi
  seed = fit_summary_params$seed
  
  # results and summary 
  sm <-  fit$summary(voi, 
                     moi,
                     qs= ~posterior::quantile2(.,probs=probsoi)) #%>% 
    #add_column(invq_real = invq( rtheta , 
     #                            fit$draws(voi, format= "draws_array"))) %>% 
    #add_column( invq_g0 = 1-invq(rep(0, length(rtheta)), 
     #                            fit$draws(voi, format="draws_array"))) 
  
  perf <- as.list(c( sum(fit$summary(c("log_lik"), mean)$mean),
                     sum(fit$summary(c("log_lik_test"), mean)$mean),
                     as.vector(t(fit$loo()$estimates))[1:4],
                     as.vector(table(cut(fit$loo()$diagnostics$pareto_k, 
                                         breaks=c(-Inf,0.5, 0.7, 1, Inf)))),
                     prmse_theta(fit, "beta", standat$beta ),
                     prmse_theta0(fit, "beta", standat$beta ),
                     prmse_thetan0(fit, "beta", standat$beta ),
                     prmse_theta_pp( fit , "beta", standat$beta),
                     prmse_theta_pp0( fit , "beta", standat$beta),
                     prmse_theta_ppn0( fit , "beta", standat$beta),
                     #prmse_theta_pp_n0(fit, "R2D2_R2", as.numeric(rtheta["R2"]) ),
                     #prmse_theta_pp(fit, "sigma", standat$sigma ),
                     prmse(standat,fit), 
                     sum(standat$beta==0),
                     sum(standat$beta!=0))
  )
  
  names(perf) <- c("lpd_train", "lpd_test",
                   "elpd_loo","elpd_loo_sd","p_loo", "p_loo_sd",
                   "paretok_good", "paretok_ok","paretok_bad", "paretok_verybad",
                   "rmse_b", "rmse_b0", "rmse_bn0",
                   "rmse_bpp", "rmse_bpp0", "rmse_bppn0",
                   #"rmse_R2",
                   #"rmse_sigma", 
                   "rmse1_train", "rmse2_train", "rmse3_train", "rmse4_train",
                   "rmse1_test", "rmse2_test", "rmse3_test", "rmse4_test",
                   "p0", "pn0")
  
  final_result <- list(seed= seed,
                       time= fit$time(), 
                       sm= sm,
                       perf= perf, 
                       rtheta= rtheta)
  
  
  
  final_result
}


