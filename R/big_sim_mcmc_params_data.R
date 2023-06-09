# Script where data corresponding to the different mcmc methods should be entered.


get_sim_mcmc_params_data <- function(R2D2_alpha_function, 
                                     gen_coef){
  
  
  # R2D2s
  R2D2_mean_R2= 0.1
  R2D2_prec_R2= 1
  api= 0.5
  
  alpha_params <- list(api=api)
  
  mcmc_params_r2d2 <- list(R2D2_mean_R2=  R2D2_mean_R2,
                           R2D2_prec_R2=  R2D2_prec_R2,
                           gen_R2D2_alpha_function = R2D2_alpha_function,
                           gen_R2D2_alpha_params = alpha_params,
                           iter=iter)
  
  #group r2d2 
  
  mcmc_params_r2d2_grouped <- list(R2D2_mean_R2=  R2D2_mean_R2,
                                   R2D2_prec_R2=  R2D2_prec_R2,
                                   gen_R2D2_alpha_function = R2D2_alpha_function,
                                   api_coefs= 0.5 , #TODO: should come from somewhere else? adaptive?
                                   api_groups= 1 ,
                                   iter=iter)
  
  # names of the mcmc_params functions
  # create the parameters given to the stan fit
  # see stan_fits.R
  
  names_list = c("r2d2_grouped")
  
  #Fits to be used
  # Stan fits that will be used
  fits_list = c(rep("r2d2grouped",1))
  
  nfits = length(fits_list) #number of fits
  
  
  mcmc_params <-  list(#mcmc_params_r2d2= mcmc_params_r2d2, 
                       mcmc_params_r2d2_grouped = mcmc_params_r2d2_grouped)
  #Used in main script of simulation
  
  fits_params <-  list(nfits=nfits, 
                       names_list = names_list,
                       nnames = length(names_list),
                       fits_list= fits_list, 
                       mcmc_params= mcmc_params)
  
  fits_params
  
}

