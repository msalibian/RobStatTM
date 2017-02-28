# test script

# R CMD INSTALL --preclean --clean robustbroli 


library(robustbroli)
library(pyinit)
data(coleman)
set.seed(123)
# source('R/lmrob2.R')

## Default for a very long time:
m1 <- lmrob2(Y ~ ., data=coleman)


