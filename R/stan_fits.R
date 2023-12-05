# TODO: Decide if sigma should be included or not? mmmm 
# at the moment it is passed as data

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
  R2D2_per_group_alphas = params$R2D2_per_group_alphas
  seed= params$seed
  
  dat <- list(
    N=n,
    p=p+1,
    X=  cbind(rep(1,n), X), 
    y= as.numeric(y), 
    G= params$groups_n, 
    pg= params$group_ps,
    Ntest= ntest, 
    Xtest= cbind(rep(1,ntest),Xtest),  
    ytest= as.numeric(ytest),  
    sigma= params$sigma, 
    R2D2_mean_R2= R2D2_mean_R2,
    R2D2_prec_R2= R2D2_prec_R2,
    R2D2_alpha_groups = R2D2_alpha_groups,
    R2D2_per_group_alphas = R2D2_per_group_alphas,
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
  pg = params$pg
  seed= params$seed
  
  dat <- list(
    N=n, 
    p=p+1,
    X=  cbind(rep(1,n), X), 
    Y= as.numeric(y), 
    G= params$groups_n,
    Ntest= ntest, 
    Xtest= cbind(rep(1,ntest),Xtest),  
    Ytest= as.numeric(ytest),  
    sigma= params$sigma, 
    bg = bg, 
    ag = ag,
    pg= pg
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
      adapt_delta=0.95), silent= FALSE)
  
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


horseshoefit <- function(params){
  
  mod <- readRDS("stan/horseshoecmdstanmodel")
  p= params$p
  n= params$n
  X= params$X
  y= params$y
  ntest= params$ntest
  ytest= params$ytest
  Xtest= params$Xtest
  
  seed= params$seed
  
  dat <- list(
    N=n, p=p+1,
    X=  cbind(rep(1,n), X), y=as.numeric(y), 
    Ntest=ntest, 
    Xtest= cbind(rep(1,ntest),Xtest),  
    ytest= as.numeric(ytest),  
    sigma= params$sigma,
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

rhorseshoefit <- function(params){
  
  mod <- readRDS("stan/horseshoecmdstanmodel")
  p= params$p
  n= params$n
  X= params$X
  y= params$y
  ntest= params$ntest
  ytest= params$ytest
  Xtest= params$Xtest
  
  seed= params$seed
  
  hs_df =  params$hs_df
  hs_df_global = params$hs_df_global
  hs_df_slab = params$hs_df_slab
  hs_scale_slab = params$hs_scale_slab
  p0 = params$p0
  
  dat <- list(
    N=n, p=p+1,
    X=  cbind(rep(1,n), X), y=as.numeric(y), 
    Ntest=ntest, 
    Xtest= cbind(rep(1,ntest),Xtest),  
    ytest= as.numeric(ytest),  
    hs_df = hs_df,
    hs_df_global = hs_df_global,
    hs_df_slab = hs_df_slab,
    hs_scale_slab = hs_scale_slab,
    p0 = p0,
    sigma= params$sigma,
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


