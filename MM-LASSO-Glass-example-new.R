

# # ## Only run once
# # # Install pense package
# library(devtools)
# install_github("dakep/pense-rpkg", ref = "develop", force=TRUE)
# #
# # # Install pyinit package
# # library(devtools)
# # install_github("dakep/pyinit", ref = "master", force=TRUE)

# data(glassVessels, package='locout')

# load('glassVessels.RData')
# x <- glass$data
# y <- glass$chemical[,'PbO']
# x <- x[, 15:500]
# library(pense)
# 
# set.seed(1)
# b0 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
# b2 <- pensem(b0, alpha=1, nlambda=30, ncores=4)
# coef(b2)[ coef(b2) != 0 ]
# which( coef(b2) != 0 )
# 
# set.seed(123)
# b0.2 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
# b2.2 <- pensem(b0.2, alpha=1, nlambda=30, ncores=4)
# coef(b2.2)[ coef(b2.2) != 0 ]
# which( coef(b2.2) != 0 )
# 
# set.seed(123456)
# b0.3 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
# b2.3 <- pensem(b0.3, alpha=1, nlambda=30, ncores=4)
# coef(b2.3)[ coef(b2.3) != 0 ]
# which( coef(b2.3) != 0 )
# 
# set.seed(17)
# b0.4 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
# b2.4 <- pensem(b0.4, alpha=1, nlambda=30, ncores=4)
# coef(b2.4)[ coef(b2.4) != 0 ]
# which( coef(b2.4) != 0 )
# 
# set.seed(29)
# b0.5 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
# b2.5 <- pensem(b0.5, alpha=1, nlambda=30, ncores=4)
# coef(b2.5)[ coef(b2.5) != 0 ]
# which( coef(b2.5) != 0 )
# 
# set.seed(654321)
# b0.6 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
# b2.6 <- pensem(b0.6, alpha=1, nlambda=30, ncores=4)
# coef(b2.6)[ coef(b2.6) != 0 ]
# which( coef(b2.6) != 0 )
# 
# set.seed(9311)
# b0.7 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
# b2.7 <- pensem(b0.7, alpha=1, nlambda=30, ncores=4)
# coef(b2.7)[ coef(b2.7) != 0 ]
# which( coef(b2.7) != 0 )
# 
# # M (10 additional) runs of the glass example
# 
# 
# load('glass-runs-1-7.RData')
# 
# M <- 10
# library(pense)
# b0s <- b2s <- vector('list', M)
# seed0 <- 123
# 
# for(j in 1:M) {
#   set.seed(seed0 + 17*j)
#   b0s[[j]] <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
#   b2s[[j]] <- pensem(b0s[[j]], alpha=1, nlambda=30, ncores=4)
#   print(c(j, which( coef(b2s[[j]]) != 0 ) ) )
#   fi <- paste('glass-runs-1-', 7+j, '.RData', sep='')
#   save.image(file=fi)
# }

# seeds

# set.seed(1)
# set.seed(123)
# set.seed(123456)
# set.seed(17)
# set.seed(29)
# set.seed(654321)
# set.seed(9311)
# 123 + 17*j, j = 1, ..., 10



# process results
b0s.b <- b0s
b2s.b <- b2s
b0s <- b2s <- vector('list', 17)
b0s[[1]] <- b0
b0s[[2]] <- b0.2
b0s[[3]] <- b0.3
b0s[[4]] <- b0.4
b0s[[5]] <- b0.5
b0s[[6]] <- b0.6
b0s[[7]] <- b0.7
b2s[[1]] <- b2
b2s[[2]] <- b2.2
b2s[[3]] <- b2.3
b2s[[4]] <- b2.4
b2s[[5]] <- b2.5
b2s[[6]] <- b2.6
b2s[[7]] <- b2.7
for(j in 1:10) {
  b0s[[j + 7]] <- b0s.b[[j]]
  b2s[[j + 7]] <- b2s.b[[j]]
}

library(pense)
for(j in 1:17) 
  print(which(coef(b2s[[j]])!=0))
  
  