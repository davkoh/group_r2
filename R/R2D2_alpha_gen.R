#--- R2D2_alpha generating functions

R2D2_alpha_api <- function(alpha_params){
  p= alpha_params$p
  api= alpha_params$api
  R2D2_alpha= rep(api, p)
  R2D2_alpha
}

R2D2_alpha_good_spikes <- function(alpha_params){
  p = alpha_params$p
  R2D2_alpha= c(rep(10,5), rep(0.5,p-10) , rep(5,2), rep(0.5,3))
  #R2D2_alpha[1: floor(0.05*p)]= 10
  #R2D2_alpha[c(1:5, p-4, p-5 )]= 10
  R2D2_alpha
}


R2D2_alpha_bad_spikes <- function(alpha_params){
  p = alpha_params$p
  #R2D2_alpha = rep(0.5, p)
  #R2D2_alpha[(floor(0.05*p)+1):p]=10
  R2D2_alpha= c(rep(.1,5), rep(5,p-10) , rep(.1,2), rep(5,3))
  R2D2_alpha
}