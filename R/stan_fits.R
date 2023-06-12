r2d2fit <- function(params){
  #file <- file.path("stan", "R2D2.stan")
  #mod <- cmdstan_model(file)
  mod <- readRDS("stan/r2d2cmdstanmodel")
  p= params$p
  n= params$n
  X= params$X
  y= params$y
  ntest= params$ntest
  ytest= params$ytest
  Xtest= params$Xtest
  
  #fit-specific params
  R2D2_mean_R2=  params$R2D2_mean_R2
  R2D2_prec_R2=  params$R2D2_prec_R2
  R2D2_alpha = params$R2D2_alpha
  seed= params$seed
  
  dat <- list(
    N=n, p=p+1,
    X=  cbind(rep(1,n), X), y=as.numeric(y), 
    Ntest=ntest, 
    Xtest= cbind(rep(1,ntest),Xtest),  ytest= as.numeric(ytest),  
    sigma= params$sigma, 
    R2D2_mean_R2= R2D2_mean_R2,
    R2D2_prec_R2= R2D2_prec_R2,
    R2D2_alpha= R2D2_alpha,
    prior_only=0
  )
  
  fit = try(
    mod$sample(
      data = dat,
      seed= seed,
      chains = 1,
      refresh = 250,
      iter_sampling =params$iter,
      output_dir= params$temp_directory,
      adapt_delta=0.99), silent= FALSE)
  
  if (is(fit, "try-error")) {
    fit.error <- TRUE
    fit <- NULL
  } else {
    fit.error <- fit$return_codes()  
  }
  
  
  if(fit.error){
    #In case there is an error
    return(list(fit.error=fit.error))
  }else{
    return(list(fit= fit, fit.error=fit.error))
  }
}


r2d2groupedfit <- function(params){
  #file <- file.path("stan", "R2D2.stan")
  #mod <- cmdstan_model(file)
  mod <- readRDS("stan/r2d2groupedcmdstanmodel")
  
  p= params$p
  X= params$X
  n= params$n
  y= params$y
  ntest= params$ntest
  ytest= params$ytest
  Xtest= params$Xtest
  
  #fit-specific params
  R2D2_mean_R2=  params$R2D2_mean_R2
  R2D2_prec_R2=  params$R2D2_prec_R2
  R2D2_alpha_groups = params$R2D2_alpha_groups
  R2D2_groups_alphas = params$R2D2_groups_alphas
  seed= params$seed
  
  dat <- list(
    N=n, p=p+1,
    X=  cbind(rep(1,n), X), 
    y= as.numeric(y), 
    G= params$groups_n, pg= params$group_ps,
    Ntest= ntest, 
    Xtest= cbind(rep(1,ntest),Xtest),  
    ytest= as.numeric(ytest),  
    sigma= params$sigma, 
    R2D2_mean_R2= R2D2_mean_R2,
    R2D2_prec_R2= R2D2_prec_R2,
    R2D2_alpha_groups = R2D2_alpha_groups,
    R2D2_groups_alphas = R2D2_groups_alphas,
    prior_only=0
  )
  
  fit = try(
    mod$sample(
      data = dat,
      seed= seed,
      chains = 1,
      refresh = 250,
      iter_sampling =params$iter,
      output_dir= params$temp_directory,
      adapt_delta=0.99), silent= FALSE)
  
  if (is(fit, "try-error")) {
    fit.error <- TRUE
    fit <- NULL
  } else {
    fit.error <- fit$return_codes()  
  }
  
  if(fit.error){
    #In case there is an error
    return(list(fit.error=fit.error))
  }else{
    return(list(fit= fit, fit.error=fit.error))
  }
}


giggfit <- function(params){
  #file <- file.path("stan", "R2D2.stan")
  #mod <- cmdstan_model(file)
  mod <- readRDS("stan/giggcmdstanmodel")
  
  p= params$p
  X= params$X
  n= params$n
  y= params$y
  ntest= params$ntest
  ytest= params$ytest
  Xtest= params$Xtest
  
  #fit-specific params
  ag =  params$ag
  bg =  params$bg
  seed= params$seed
  
  dat <- list(
    N=n, 
    p=p+1,
    X=  cbind(rep(1,n), X), 
    Y= as.numeric(y), 
    G= params$groups_n,
    pg= params$pg,
    Ntest= ntest, 
    Xtest= cbind(rep(1,ntest),Xtest),  
    Ytest= as.numeric(ytest),  
    sigma= params$sigma, 
    bg = bg, 
    ag = ag
    #prior_only=0
  )
  
  fit = try(
    mod$sample(
      data = dat,
      seed= seed,
      chains = 1,
      refresh = 250,
      iter_sampling =params$iter,
      output_dir= params$temp_directory,
      adapt_delta=0.99), silent= FALSE)
  
  if (is(fit, "try-error")) {
    fit.error <- TRUE
    fit <- NULL
  } else {
    fit.error <- fit$return_codes()  
  }
  
  if(fit.error){
    #In case there is an error
    return(list(fit.error=fit.error))
  }else{
    return(list(fit= fit, fit.error=fit.error))
  }
}


