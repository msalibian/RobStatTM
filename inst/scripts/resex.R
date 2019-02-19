# resex.R
# Example 8.6
# Figures 8.12, 8.13
# Table 8.5 

# Must install robustarima
library(robustarima)

# resex=scan("resex0.txt") uses wrong .txt file, correct code below
# resex=scan("resex.txt") replaced below
data(resex,package='RobStatTM')
 
#filered tau
resxar=arima.rob(formula = resex~ 1, p = 2, sd = 1, sfreq = 12)
 
#autoregressive coefficients
arcoeftau=resxar$model$ar
arcoeftau
#mean of the the differenced series
meantau=resxar$regcoef
names(meantau)="mean"
meantau
#intercept of the differenced series
# intercepttau=meantau*(1-sum(arcoef)) incorrect arcoef, correct code below
intercepttau=meantau*(1-sum(arcoeftau))
names(intercepttau)="intercept"
# intercept # no object intercept, correct code below
intercepttau
#sorted innovations
innovtau=sort(abs(resxar$innov[15:89]))

#Figure 8.22 
  
plot(1:89, resex, type="l",xlab="index", ylab="RESEX")
points(1:89, resxar$y.robust )
 
#LS
 
sresx=resex[13:89]-resex[1:77]
resxls=lm(formula = sresx[3:77] ~ sresx[2:76] + sresx[1:75])
#autoregressive coefficients
arcoefls=resxls$coef[2:3]
names(arcoefls)=c("AR(1)", "AR(2)")
arcoefls
#intercept
interceptls=resxls$coeff[1]
names(interceptls)="intercept"
interceptls

#sorted innovations
innovls=sort(abs(resxls$residuals))

# Figure 8.23 
    
ttt=ppoints(72)
plot(ttt,innovls[1:72], type="l", xlab="probability",ylab="quantiles")
lines(ttt,innovtau[1:72],lty=1)
#lines(ttt,resxgminnov[1:72],lty=1)
 
#text(c(.5,.5,.5),c(.83,1,1.25),c("TAU","GM","LS"))
text(c(.5,.5),c(.83,1.25),c("TAU","LS"))

 
 
 
