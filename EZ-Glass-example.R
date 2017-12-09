# Compare with Ezequiel's MM-LASSO
library(mmlasso)
load('glassVessels.RData')
x <- glass$data
# table(glass$group)
# matplot(t(x), col='gray', type='l', lty=1)
y <- glass$chemical[,'PbO']
x <- x[, 15:500]

M <- 17
mms <- vector('list', M)
for(j in 1:M) {
  set.seed(123 + 17*j)
  mms[[j]] <- mmlasso(x=x, y=y, ncores=4)
  print(c(j, which(mms[[j]]$coef.MMLasso!=0)))
  print(c(j, which(mms[[j]]$coef.MMLasso.ad!=0)))
}

save.image('glass-ez-runs.RData')
          
