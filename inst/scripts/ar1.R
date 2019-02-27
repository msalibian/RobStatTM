# ar1.R
# Simulated AR1  with AO and IO
# AR(1) with AO and IO
# Figure 8.6

# Must install robustarima
library(robustarima)

set.seed(1000)
n.innov = 200
n= 100
 phi=0.9

n.start = n.innov - n
innov = rnorm(n.innov)

x= arima.sim(model = list(ar = phi), n, innov = innov, n.start = n.start)

#10% of additive outliers of equistant  outliers

ao=rep(0,n)
tt=seq(10,100,10)
ao[tt]=4
 xAO=x +ao
#innovation outlier at observation
xIO=x
xIO[1:49]=x[1:49]
xIO[50]=phi*xIO[49]+10
u=rnorm(50)
for (i in 51:100)
xIO[i]=phi*xIO[i-1] +u[i-50]

par(mfrow=c(3,1))
plot(x)
title("Gaussian AR(1) series without outliers")
plot(xAO)
points(tt,xAO[tt],pch=1)
title("Gaussian AR(1) series with 10% of additive outliers")
plot(xIO)
points(50,xIO[50] ,pch=1 )
title("Gaussian AR(1) series with one innovation  outlier")
par(mfrow=c(1,1))







