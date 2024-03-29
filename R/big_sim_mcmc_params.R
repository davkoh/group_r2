#MCMC functions

#---- Make the MCMC params

make_mcmc_params <- function(name, mcmc_params, data_gen_params){
  get(paste0(name, "_mcmc_params"))(mcmc_params, data_gen_params)
}



#---- MCMC params

#--- R2D2
r2d2_mcmc_params <- function(mcmc_params, data_gen_params){
  #fit-specific params
  R2D2_mean_R2=  mcmc_params$R2D2_mean_R2
  R2D2_prec_R2=  mcmc_params$R2D2_prec_R2
  
  alpha_params= mcmc_params$gen_R2D2_alpha_params
  alpha_params$p=data_gen_params$p
  
  gen_R2D2_alpha <- get( mcmc_params$gen_R2D2_alpha_function, mode = "function")
  
  R2D2_alpha = gen_R2D2_alpha(alpha_params)
  
  return(list( R2D2_mean_R2= R2D2_mean_R2,  R2D2_prec_R2= R2D2_prec_R2 , R2D2_alpha= R2D2_alpha  ))
  
}

#--- grouped R2D2

r2d2_grouped_mcmc_params <- function(mcmc_params, data_gen_params){
  #fit-specific params
  R2D2_mean_R2=  mcmc_params$R2D2_mean_R2
  R2D2_prec_R2=  mcmc_params$R2D2_prec_R2
  
  alpha_params= mcmc_params$gen_R2D2_alpha_params
  alpha_params$G= data_gen_params$groups_n
  
  gen_R2D2_alpha <- get( mcmc_params$gen_R2D2_alpha_function, 
                         mode = "function")
  
  
  # First decomposition. 
  # group alpha
  
  alpha_params$p = data_gen_params$groups_n
  alpha_params$api= mcmc_params$api_groups
  R2D2_alpha_groups = gen_R2D2_alpha(alpha_params)
  
  
  # Second decomposition
  # alpha inside groups
  alpha_params$p= data_gen_params$p
  alpha_params$api= mcmc_params$api_coefs
  R2D2_per_group_alphas = gen_R2D2_alpha(alpha_params)
  
  return(list( R2D2_mean_R2= R2D2_mean_R2, 
               R2D2_prec_R2= R2D2_prec_R2 , 
               #R2D2_alpha= R2D2_alpha ,
               R2D2_alpha_groups = R2D2_alpha_groups,
               R2D2_per_group_alphas = R2D2_per_group_alphas))
  
}

#--- gigg prior

gigg_mcmc_params <- function(mcmc_params= NULL, data_gen_params){
  #fit-specific params
  groups_n <- data_gen_params$groups_n
  p <- data_gen_params$p
  #TODO: Discuss with DK
  
  # Horsehoe gigg see paper

  ag <- rep(0.5, groups_n )
  bg <- rep(0.5, p)
  pg <- rep( c(1:(groups_n)), each= p/groups_n )

  return(list( ag= ag, 
               bg= bg, 
               pg= pg))
  
}


#---- Horseshoe

horseshoe_mcmc_params <- function(mcmc_params, data_params){
  #fit-specific params
  
  return(list( dummy= NULL ))
  
}


#--- Regularized Horseshoe


rhorseshoe_mcmc_params <- function(mcmc_params, data_params){
  #fit-specific params
  
  hs_df =  mcmc_params$hs_df
  hs_df_global = mcmc_params$hs_df_global
  hs_df_slab = mcmc_params$hs_df_slab
  hs_scale_slab = mcmc_params$hs_scale_slab
  p0 = data_params$p*0.5
  
  
  return(list(hs_df = hs_df, 
              hs_df_global= hs_df_global, 
              hs_df_slab= hs_df_slab, 
              hs_scale_slab= hs_scale_slab, 
              p0= p0))
  
}

