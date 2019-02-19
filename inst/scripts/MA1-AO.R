# MA1-AO.R
# Example 8.5
# Figure 8.11, Table 8.4
# Robust fitting of a simulated MA(1) series

# Must install robustarima	
library(robustarima)

set.seed(200)
n.innov = 300
n = 200
theta=-0.8

n.start = n.innov - n 
innov = rnorm(n.innov)
n.start = n.innov - n 
 
ma1 = arima.sim(model = list(ma = theta), n=n, innov = innov, n.start = n.start)

#ma1 <- arima.sim(model=list(ma=c(-.8)),n=200,n.innov=n.innov 
mac=ma1
mac[20*(1:10)]=ma1[20*(1:10)]+4
ma1tau=arima.rob(mac~1,q=1)
ma1tau
ma1ls=arima(mac,order=c(0,0,1),method="CSS")

#Figure 8.21
 
plot(1:200, mac, ylab='series', xlab='index', type='l', lty=1 )
outma=seq(20,200,20)
points(outma, mac[outma])
lines(1:200, ma1tau$y.robust, lty=2 )
 
 



 


 
