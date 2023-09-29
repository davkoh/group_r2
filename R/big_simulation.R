# Big simulation test 

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
library(mvtnorm)
library(MASS)
library(ParallelLogger)

# NOTE: Additionally required packages: 
# library(posterior)
# library(sirt)
# library(expm)

# not actually required packages
# library(tidyverse)
# library(foreach)
# library(iterators)
# library(compositions)
# library(ggplot2)
# library(corrplot)
# library(ggcorrplot)
# library(bayesplot)
# library(paletteer)# Simulation Design

#--- Functions

source("R/aux_functions.R")
source("R/stan_fits.R")
source("R/big_sim_full_simulation.R")
source("R/big_sim_mcmc_params.R")
source("R/big_sim_mcmc_params_data.R")
source("R/R2D2_alpha_gen.R")


#--- 

n= 500 #size of training data
ntest= 200 # size of test data 
ps= c(50) #num overall of coefficients
rho_in = c(0.8) #correlation of X, X ~ MVN
rho_out = c(0.2)
nus= c(0.75) #sparsity level 
type="AR" #type of correlation matrix used
alpha= 0 #real intercept

Gs= c(5) #number of groups

R2= c(0.1,0.3, 0.5,0.7)  # signal to noise ratio
type <- c("dist", "con")

#simulation conditions
sim_conds <- expand.grid(n=n, ntest=ntest, p=ps, 
                         rho_in= rho_in, rho_out= rho_out, nu=nus, type=type, 
                         alpha= alpha, Gs= Gs, R2 = R2)

sim_conds <-  sim_conds %>% 
              add_column( id= seq_len(nrow(sim_conds)), .before=1)


#--- Extra conditions

#Function to generate concentration vector 
R2D2_alpha_fns = c("R2D2_alpha_api")   

gen_coef_fns= c( "gen_coef_sparse_norm")

gen_coef_params <-  NULL

extra_conds <- data.frame( R2D2_alpha_function= R2D2_alpha_fns, 
                           gen_coef= gen_coef_fns )


#--- Simulation specific params

# TODO: change this to the desired number of simulations

nsims = 50 # number of simulation per condition
iter= 2000 # MCMC iters

# sm: summary quantities of interest
# qoi: quantities of interest
# voi: variables of interest
# moi: metrics of interest
# probsoi: probabilities of interest

smqoi <- list(voi= c("beta",
                     "b_Intercept"), 
                     #"R2",
                     #"sigma"), 
                     #"R2D2_phi", #
                     #"lambdas") ,
              moi=c("mean",
                    "sd", 
                    "median", 
                    "mad",
                    "ess_basic",
                    "ess_bulk",
                    "ess_tail",
                    "rhat_basic", 
                    "rhat"),
              probsoi= c(0.01, 0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.33, 0.40, 0.5, 0.6,  0.67,
                         0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99))

#--- Run simulation

path_results <- "R/big_simulation/final_results"
global_seed =  117

# TODO: change to your cmdstan path
my_cmd_stan_path = cmdstanr::cmdstan_path()
cmdstanr::set_cmdstan_path(my_cmd_stan_path)

ncores_simulation = parallel::detectCores()

temp_stan_directory <- paste0("temp_stan_files/")
dir.create(temp_stan_directory)

tot= nrow(extra_conds) # For each extra condition we run all of the conditions of the simulations
my_name= "new_GGG" #Name of the files. Format my_name_global_seed_sim_cond_id 

#----- Compile fits that will be used


file <- file.path("stan", "R2D2.stan")
mod_R2D2<- cmdstan_model(file)
saveRDS(mod_R2D2, "stan/r2d2cmdstanmodel") 

file <- file.path("stan", "R2D2_grouped.stan")
mod_R2D2_grouped <- cmdstan_model(file)
saveRDS(mod_R2D2_grouped, "stan/r2d2groupedcmdstanmodel") 

file <- file.path("stan", "horseshoe.stan")
mod_horseshoe<- cmdstan_model(file)
saveRDS(mod_horseshoe, "stan/horseshoecmdstanmodel")


file <- file.path("stan", "rhorseshoe.stan")
mod_rhorseshoe<- cmdstan_model(file)
saveRDS(mod_rhorseshoe, "stan/rhorseshoecmdstanmodel")

file <- file.path("stan", "gigg.stan")
mod_gigg <- cmdstan_model(file)
saveRDS(mod_gigg, "stan/giggcmdstanmodel") 

stan_models_list <- list(mod_R2D2= mod_R2D2,
                         mod_R2D2_grouped= mod_R2D2_grouped,
                         mod_horseshoe= mod_horseshoe,
                         mod_rhorseshoe = mod_rhorseshoe,
                         mod_gigg = mod_gigg)

#- Run simulation


for(i in 1:tot){
 
  #TODO: CHANGE SIMULATION so that conditions can vary per group. 
  
  # Conditions to be given to other script
  gen_coef =  extra_conds$gen_coef[i] # how to generate coefficients inside all groups
  
  R2D2_alpha_function= extra_conds$R2D2_alpha_function[i] #how is R2D2 alpha generated?
  exp_name= paste0(my_name, as.character(i),"_", global_seed) #Experiment name
  
  # This is a function
  # TODO: eventually, someday with some energy, change this and make efficient
  fits_params <-  get_sim_mcmc_params_data(R2D2_alpha_function,
                                           gen_coef)
  
  # simulation parameters: 
  # gencoef_function: function that generates b
  # gencoef_params: params needed to generate b
  # fits_params: params of the fits. See extra script sim_mcmc_params_data
  # nsims: number of simulations per condition
  
  sim_params <-  list(fits_params = fits_params,
                      gen_coef =extra_conds$gen_coef[i],
                      nsims= nsims, 
                      cmdstan_path= my_cmd_stan_path)
  
  # Calculate sigma based on R2, simulation conditions and extra conditions
  # sigma is a function of R2
  
  #sigma=c()
  # for(j in 1:nrow(sim_conds)){
  #   
  #   groups_n = sim_conds[j,]$Gs # number of groups
  #   group_rhos= rep(sim_conds[j,]$rho, groups_n) #rho in each group
  #   group_cor_type= rep("AR", groups_n) #type of correlation structure 
  #   #number of covariates inside each group
  #   group_ps = rep(sim_conds[j,]$p/groups_n, groups_n) 
  #   
  #   # how to generate coefficients inside each group
  #   group_nus = rep(sim_conds[j,]$nu ,groups_n)
  #   group_gen_coef_functions= rep(extra_conds$gen_coef[i], groups_n) 
  #   
  #   
  #   group_params <- list(seed= global_seed,
  #                        groups_n= groups_n,
  #                        group_nus= group_nus,
  #                        group_rhos= group_rhos ,
  #                        group_cor_type= group_cor_type,
  #                        group_ps= group_ps, 
  #                        group_gen_coef_functions= group_gen_coef_functions)
  #   
  #   
  #   # Covariance matrix of X
  #   covx <- get_sigmaX_grouped(group_params)
  #   
  #   # Find sigma given R2
  #   
  #   sigma_x <- diag(covx)
  #   sigma_b <-  gen_coef_groups_var(group_params)
  #   sigma[j] <- sqrt(sum(sigma_x*sigma_b)*(1-sim_conds[j,]$R2)/sim_conds[j,]$R2)
  #   
  # }
  # 
  # sim_conds$sigma = sigma
  # 

  # RUN FULL SIMULATION! 
  
  full_simulation(sim_conds= sim_conds,  
                  sim_params= sim_params,
                  smqoi= smqoi,
                  ncores_simulation = ncores_simulation,
                  seed= global_seed, 
                  path = path_results, 
                  special_name= exp_name )
  
  
  
}


#--- 

# 
# temp <-   full_simulation(sim_conds= sim_conds,  
#                           sim_params= sim_params,
#                           smqoi= smqoi,
#                           ncores_simulation = ncores_simulation,
#                           seed= global_seed, 
#                           path = path_results, 
#                           special_name= exp_name )
# 
# sim_cond <- sim_conds[1,]
# 
# temp <- cond_sim(sim_params,
#          sim_cond, 
#          smqoi, seed = 1)
# 
# temp
# 


# how to make better groupings?

# If we have prior knowledge about the coeficients could we possibly add it to the prior?
# This are groups with high signal. How high should the signal be?
# Grouping via signal strength could be categorical. Priors could be the same (they shouldnt)
# Hypothesis: grouping of a signal strength could be important
# g groups of different signal strength

# Priors are noninformative but the grouping might be

# Have a whole bunch of groups are correlated and randomly assign the coefficients to the groups
# then make them large or small (decouple selection and assignation!)

# 







