dgp_boss <- function(N = 500,K = 50,G = 5,rho = 0.8,type = "dist", SNR = 0.7){
# "dist" = distributed, "conc" = concentrated
  # Generate X
  mux <- array(0,c(K,1))
  rho_in <- rho
  rho_out <- 0.2
  block <- diag(pg)+rho_in -diag(pg)*rho_in
  
  covx <- kronecker(diag(G),block)
  covx[covx==0] <- 0.2
  X <- mvrnorm(N,mux,covx)
  
  # Generate beta
  if (type == "dist"){
    beta <- t(cbind(t(rep(0.5,5)),t(rep(1,5)),t(rep(0,K-K/G))))
  } else {
    beta <- array(0,c(K,1))
    beta[seq(1,K,K/G)] <- c(0.5,1,1.5,2,2)
  }
  
  # Generate sigma
  sigma <- sqrt((1-SNR)/SNR * t(beta)%*%t(X)%*%X%*%beta / N)
  # Generate Y
  Y <- X%*%beta + rnorm(N)*as.numeric(sigma)
  
  return(list(Y = Y, X = X, beta = beta, sigma = sigma , empR2 = t(beta)%*%t(X)%*%X%*%beta / (t(beta)%*%t(X)%*%X%*%beta + N*sigma^2) ))
}
