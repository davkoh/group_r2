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
    params_temp= list(p=p, nu=nu, seed= seed)
    
    #calculate variance per group
    gen_coef_fn_var <- get( paste0(params$group_gen_coef_functions[i], "_var"),
                            mode = "function")
    #return vector of variances
    var_b=c(var_b,gen_coef_fn_var(params_temp))
    
  }
  var_b
}



