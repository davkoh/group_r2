library(gigg)
library(cmdstanr)
rm(list = ls())
X = concentrated$X
C = concentrated$C
Y = as.vector(concentrated$Y)
grp_idx = concentrated$grps

gf = gigg(X,array(1,c(dim(X)[1],1)), Y, method = "fixed", grp_idx, n_burn_in = 500, n_samples = 1000, 
          n_thin = 1, verbose = TRUE, btrick = FALSE, stable_solve = TRUE)

mod <- cmdstan_model("/Users/dk/Documents/group_r2/Stan/groupr2_dk_fixedg_v2.stan")

dat<- list( N = dim(X)[1],
            p = dim(X)[2],
            Y = Y,
            X = X,
            G = length(unique(grp_idx)),
            pg = 10,#histc(grp_idx,unique(grp_idx))[["cnt"]],
            mean_R2 = 0.33,
            prec_R2 = 3,
            cons = rep(1,length(unique(grp_idx))),
            cons_g = rep(1,dim(X)[2]),
            sigma_sd = std(Y),
            alpha_mean = 0,
            alpha_sd = 5
)


fit <- mod$sample(
  data = dat,
  seed = 123, 
  chains = 4,
  adapt_delta = 0.999,
  max_treedepth = 15,
  parallel_chains = 4,
  refresh = 100,
  iter_warmup = 2000,
  iter_sampling = 2000
)


# Quick checks
fit$summary(variables = c("psi"))
hist(fit$draws(variables = "postr22"),xlim = c(0,1))
sqrt(var(colMeans(fit$draws(variables = c("ypred"),format = "matrix"))-Y))
sqrt(var(X%*%gf$beta.hat-Y))
log_liks <- fit$draws(variables = "logliks")
loo(log_liks)

R2_gigg <- array(0,c(gf$draws$n_samples,1))

for (i in 1:gf$draws$n_samples){
  R2_gigg[i] <- t(gf$draws$betas[,i])%*%t(X)%*%X%*%gf$draws$betas[,i] /(t(gf$draws$betas[,i])%*%t(X)%*%X%*%gf$draws$betas[,i] + dim(X)[1]*gf$draws$sigma_sqs[i])
}
hist(R2_gigg,xlim = c(0,1))