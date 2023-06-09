# Small script to simulate data 

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

library(bayesplot)
library(posterior)
library(ggplot2)

#--- Functions

source("R/aux_functions.R")
source("R/stan_fits.R")

#--- Prepare data

seed=1
n= 100 #train size
ntest= 500 # test size


G= 5
groups_n = G # number of groups
group_rhos= c(0.1, -0.95, 0.5, 0.5, 0.75) #rho in each group
group_cor_type= rep("AR", G) #type of correlation structure 
group_ps = rbinom(G, 10, 0.5) #number of covariates inside each group
group_ps = rep(5, G)
p= sum(group_ps)

# how to generate coefficients inside each group
group_nus = rep(0.5,G)
group_gen_coef_functions= rep("gen_coef_sparse_norm", G) 


# list of parameters that will be used in data generating procedure
group_params <- list(seed= seed,
                     groups_n= groups_n,
                     group_nus= group_nus,
                     group_rhos= group_rhos ,
                     group_cor_type= group_cor_type,
                     group_ps= group_ps, 
                     group_gen_coef_functions= group_gen_coef_functions)


# Covariance matrix of X
covx <- get_sigmaX_grouped(group_params)

# Plot correlation matrix of X
corrplot::corrplot(cov2cor(get_sigmaX_grouped(group_params)))

# Find sigma given R2

R2= 0.5 # signal to noise ratio
nus = group_nus

sigma_x <- diag(covx)
sigma_b <-  gen_coef_groups_var(group_params)
sigma <- sqrt(sum(sigma_x*sigma_b)*(1-R2)/R2)


# Design matrix X

set.seed(seed)
X <- rmvnorm(n, mean= rep(0,p), sigma=covx)
X <- scale(X)

set.seed(seed+12313451)
Xtest <- rmvnorm(ntest, mean= rep(0,p), sigma=covx)
Xtest <- scale(Xtest)

set.seed(seed)

# coefficients
alpha <- 0
beta <- gen_coef_groups(group_params)

#vector of real parameters
real_params <- list(beta= beta, 
                    sigma= sigma, 
                    R2=R2, 
                    tau2= R2/(1-R2))

# Generate y
y= as.numeric(cbind(rep(1,n), X)%*%c(alpha,beta)+rnorm(n,0, sigma))          
ytest= as.numeric(cbind(rep(1,ntest), Xtest)%*%c(alpha,beta)+rnorm(ntest,0, sigma))          

#--- stan

#---  R2D2 grouped JA

file <- file.path("stan", "R2D2_grouped.stan")
mod_R2D2grouped <- cmdstan_model(file)
saveRDS(mod_R2D2grouped, "stan/r2d2groupedcmdstanmodel") 

#hyperparametes
R2D2_mean_R2 = 0.5
R2D2_prec_R2 =  1
R2D2_alpha_groups= rep(5, G)
R2D2_groups_alphas= rep(0.5, p)

params <- list(seed= seed, 
               p=p, n=n, X=X, y=y, 
               ntest= ntest, Xtest= Xtest, ytest= ytest,
               G= G,
               pg= group_ps,
               sigma= sigma,
               R2D2_mean_R2= R2D2_mean_R2, 
               R2D2_prec_R2= R2D2_prec_R2, 
               R2D2_alpha_groups= R2D2_alpha_groups, 
               R2D2_groups_alphas= R2D2_groups_alphas)

fit <- r2d2groupedfit(params)
fit <- fit$fit

draws_df_JA <- as_draws_df(fit)


#--- 

file <- file.path("stan", "groupr2_dk_fixedg_v2.stan")
mod_R2D2grouped_dkv2 <- cmdstan_model(file)
saveRDS(mod_R2D2grouped, "stan/r2d2groupedcmdstanmodel_dkv2") 

#hyperparametes
R2D2_mean_R2= 0.5
R2D2_prec_R2 =  1
R2D2_alpha_groups= rep(5, G)
R2D2_groups_alphas= rep(0.5, p)

dat_dk <- list(seed= seed, 
               p=p, N=n, X=X, Y=y, 
               ntest= ntest, Xtest= Xtest, ytest= ytest,
               G= G,
               pg= group_ps[1],
               sigma= sigma,
               mean_R2= R2D2_mean_R2, 
               prec_R2= R2D2_prec_R2, 
               cons= R2D2_alpha_groups, 
               cons_g = R2D2_groups_alphas,
               sigma_sd= 3,
               alpha_mean=0, 
               alpha_sd= 5)

fit_dkv2 <- mod_R2D2grouped_dkv2$sample(
  data = dat_dk,
  seed= seed,
  chains = 1,
  refresh = 250,
  adapt_delta=0.99)

draws_df_dk <- as_draws_df(fit_dkv2)

bayesplot::mcmc_recover_hist( x= fit$draws("R2D2_R2"), true= real_params$R2)+ ggtitle("JA")
bayesplot::mcmc_recover_hist( x= fit_dkv2$draws("R2"), true= real_params$R2)+ ggtitle("DK")



bayesplot::mcmc_recover_hist( x= fit_dkv2$draws("sigma"), true= real_params$sigma)


bayesplot::mcmc_recover_intervals(x=fit$draws("beta"), true=as.numeric(real_params$beta), 
                                  size=2, point_est="median")+ xaxis_text(on = FALSE)+ 
            ggtitle("JA")
              

bayesplot::mcmc_recover_intervals(x=fit_dkv2$draws("beta"), true=as.numeric(real_params$beta), 
                                  size=2, point_est="median")+ xaxis_text(on = FALSE)+
  ggtitle("DK")
    

ytilde = fit$draws("y_tilde")
ytilde = as.matrix(as.data.frame(ytilde))

ytilde_dkv2 = fit_dkv2$draws("ypred")
ytilde_dkv2 = as.matrix(as.data.frame(ytilde_dkv2))

I=  sample(1:  nrow(ytilde) , 100,  replace=FALSE) 


ppc_dens_overlay(y, ytilde_dkv2[I,])+ggtitle("DK")
ppc_dens_overlay(y, ytilde[I,])+ggtitle("JA")



lambdas <- fit$draws("lambdas") 
mcmc_areas_ridges(lambdas)+ggtitle("JA")


print(fit$summary(variables="beta"), n=30)
print(fit_dkv2$summary(variables="beta"), n=30)
real_params$beta

