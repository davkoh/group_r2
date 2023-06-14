#notes

# 1. Nice to pregenerate the data  (TODO)
# gps --- put in a list 
# seed might change and data might be slightly different 
# you want to make sure that the data is the same


# good practice: pregenerate the data
# script - loop - gen data - load the data


# Q: id_nsim_seed unique identifier? 

# 2. Implementations 

# beta, sigma and r2
# Vary sigma and have control over R2
# r2 and sigma and then infer sigma

# sigma, beta -- r2. p grows r2 

# r2, beta --- sigma --- posterior inference and lm( y ~ x)
# conditional r2 != population r2 bc n is finite
# since n is finite then

# conditional r2, beta -- sigma per simulation per condition (TODO!)

# How to solve for beta? 
# 1. same across the vector - we need some sparsity? 
# 2. maybe looking at the original r2d2 paper?

# Start and imitate what they do. 
# In the end we will compare with them ---- 


# r2d2 - gigg 
# 

# 3. Change sparsity 

# Simplex on the groups -- fixed and known - DGP reflects this
# Only second simplex moves


# 4. Interesting to have the energy score, 
# Take joint cdf and evaluate the points of each draw 
# rmse
# david will share the energy score

# 5. 
# Are we shrinking group wise correctly? Gains or losses 
#





