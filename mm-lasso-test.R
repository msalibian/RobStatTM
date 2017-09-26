library(devtools)
install_github("dakep/pense-rpkg", ref = "develop")


library(pense)
library(mmlasso)
n <- 80
p <- 50
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 5), rep(0, p-5))) + rnorm(n, sd=.5)

a <- mmlasso(x=x, y=y, ncores=2)
b1 <- pense(X=x, y=y, alpha=0, nlambda=30, ncores=2)
b2 <- pensem(b1, alpha=1, nlambda=30, ncores=2)
b3 <- pensem(b1, alpha=0, nlambda=30, ncores=2)
a$coef.MMLasso
as.vector(b1$coefficients[, b1$lambda==b1$lambda_opt])
as.vector(b2$coefficients[, b2$lambda==b2$lambda_opt])
b2$lambda
as.vector(b3$coefficients[, b3$lambda==b3$lambda_opt])

b3 <- pensem(b1, alpha=0, nlambda=30, ncores=2)
# Error in checkForRemoteErrors(val) : 
#   one node produced an error: BLAS/LAPACK routine 'DLASCL' gave error code -4
b4 <- pensem(x=x, y=y, alpha=0, nlambda=30, ncores=2)
# Error in checkForRemoteErrors(val) : 
#   one node produced an error: BLAS/LAPACK routine 'DLASCL' gave error code -4
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# 48: In pen_s_reg(std_data$xs, std_data$yc, alpha = alpha,  ... :
# PENSE did not converge for lambda = 7.845403
# 49: In pen_s_reg(std_data$xs, std_data$yc, alpha = alpha,  ... :
# PENSE did not converge for lambda = 4.872128
                                  

as.vector(b4$coefficients[, b4$lambda==b4$lambda_opt])



a$coef.MMLasso
as.vector(b2$coefficients[, b2$lambda==b2$lambda_opt])
b2$lambda
as.vector(b3$coefficients[, b3$lambda==b3$lambda_opt])
as.vector(b4$coefficients[, b4$lambda==b4$lambda_opt])

# Fixed code
# > as.vector(b2$coefficients[, b2$lambda==b2$lambda_opt])
# [1]  0.67548491  6.94300461  0.04585328  2.56095179  0.00000000  3.89950818 -1.39608100 -1.11179883  0.00000000
# [10]  0.00000000  0.72753827  0.00000000  0.01223817  0.00000000  0.00000000 -0.39182042 -0.04672923  0.03655892
# [19] -0.06892013  0.03508527  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [28]  0.58541301 -0.28254020  0.00000000  0.00000000  1.43524046 -0.13342441  0.00000000  0.91048133  0.00000000
# [37] -0.20838049  0.00000000  0.00000000  0.00000000 -2.74201252  0.00000000  0.00000000  0.00000000  0.00000000
# [46] -0.70331155  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000

# > a$coef.MMLasso
# [1] -0.034880450  7.029282568  6.944518729  6.826255640  6.865985596  6.980971447  0.000000000  0.000000000
# [9]  0.000000000  0.000000000  0.000000000 -0.009214762  0.000000000  0.000000000  0.000000000  0.000000000
# [17]  0.000000000  0.000000000 -0.048080790  0.000000000  0.000000000  0.000000000  0.000000000 -0.029418986
# [25]  0.000000000  0.000000000 -0.039625478  0.000000000  0.000000000  0.000000000  0.000000000 -0.029157987
# [33]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
# [41]  0.091785702  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000 -0.033593114
# [49]  0.000000000  0.000000000  0.000000000
# > as.vector(b2$coefficients[, b2$lambda==b2$lambda_opt])
# [1]  8.069141e-01  7.087097e-17  1.670020e-17 -7.142006e-17  7.669765e-17 -1.046849e-16 -2.424006e-17
# [8]  7.822941e-18 -8.158476e-17  4.360117e-17  1.825734e-16  5.608968e-18  5.908426e-18  3.595310e-18
# [15] -2.834246e-17  6.885836e-17 -4.495084e-17  4.896665e-17  5.839672e-17 -2.801863e-17 -3.088070e-17
# [22] -4.299208e-17 -3.767627e-17 -2.309696e-17 -4.352130e-17 -6.031491e-17 -4.609111e-17  9.072249e-17
# [29] -9.053639e-17 -3.928071e-17 -1.008295e-16 -2.040313e-17 -3.953020e-18  1.775933e-17 -1.472012e-17
# [36]  3.226361e-17 -8.631586e-17 -2.726083e-17  7.206761e-17  1.664325e-17  5.651056e-17 -2.383603e-17
# [43] -4.079189e-17 -6.232289e-17 -3.316202e-17 -4.742614e-18 -9.502352e-17  1.655001e-17 -1.297285e-16
# [50] -9.299888e-17  2.724369e-19
# > b2$lambda
# [1] 1.628094e-04 2.402642e+11 3.573575e+11 5.315165e+11 7.905522e+11 1.175829e+12 1.748872e+12 2.601189e+12
# [9] 3.868884e+12 5.754393e+12 8.558808e+12 1.272996e+13 1.893393e+13 2.816143e+13 4.188596e+13 6.229917e+13
# [17] 9.266080e+13 1.378192e+14 2.049857e+14 3.048860e+14 4.534728e+14 6.744737e+14 1.003180e+15 1.492082e+15
# [25] 2.219251e+15 3.300807e+15 4.909463e+15 7.302100e+15 1.086079e+16 1.615382e+16 2.402642e+16
# > 
# 
# 
