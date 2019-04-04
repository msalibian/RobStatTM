library(RobStatTM)
options(digits=4)

#Trimmed mean
trimean<-function(x,alfa) {
  n=length(x); m=floor(n*alfa)
  xs=sort(x); mu=mean(xs[(m+1):(n-m)])
  xs=xs-mu; A=m*xs[m]^2 +m*xs[n-m+1]^2 +sum(xs[(m+1):(n-m)]^2)
  mu.std=A/(n-2*m); mu.std=sqrt(mu.std/n)
  return(list(mu=mu, mu.std=mu.std))
}

n=24; qn=qnorm(0.975)

data(flour)
x = as.vector(flour[,1])
resu = locScaleM(x,eff = 0.95)

muM=resu$mu; muMst=resu$std.mu;  h=muMst*qn
interM=c(muM-h, muM+h)

xbar=mean(x); smed=sd(x)/sqrt(n); h=smed*qn
intermean=c(xbar-h,xbar+h)

resu=trimean(x,0.25)
mu25=resu$mu; ss25=resu$mu.std; h=ss25*qn
inter25=c(mu25-h,mu25+h)

# Table 2.4
print("Mean, bisquare M- estimator, and 25% trimmed mean")
print(c(xbar,muM,mu25))
print("Their estimated standard deviations")
print(c(smed, muMst, ss25))

print("Their 0.95 confidence intervals")
print(rbind(intermean, interM, inter25))
