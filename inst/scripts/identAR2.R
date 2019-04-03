# identAR2.R
# Example 8.3
# Figure 8.7, 8.8

# Must install robustarima
library(robustarima)

set.seed(700)
n.innov = 300
n = 200
phi=c(4/3, -5/6  )

n.start = n.innov - n
innov = rnorm(n.innov)

x= arima.sim(model = list(ar = phi), n, innov = innov, n.start = n.start)
ao = ifelse(runif(n)>.1, 0, rnorm(n,4,1))
ao = sign(runif(n,-1,1))*ao
y = x + ao

no=sum(ao!=0) # Number of additive outliers in y

#Figure 8.7
par(mfrow=c(2,1))
plot(x,  ylab=expression(x[t]),ylim=c(-9,9))
plot(y,  ylab=expression(y[t]),ylim=c(-9,9))
ao.times = (1:n)[ao != 0]
points(ao.times, y[ao != 0])
par(mfrow=c(1,1))

# Robust automatic AR(p) fit
out=arima.rob(y~1, auto.ar=TRUE)

# Warning message, that optimization in arima.rob did not converge,
# did not effect the following code results

#Figure 8.8
par(mfrow=c(4,2))
acf1=acf(x ,10, plot=FALSE )
plot(acf1, main="acf of x")
acf2=acf(x,10,"partial", plot=FALSE)
plot(acf2, main="pcf of x")
acf3=acf(y,10, plot=FALSE)
plot(acf3, main="acf of y")
acf4=acf(y,10,"partial", plot=FALSE)
plot(acf4, main="pcf of y")

# Procedure (A), out$y.robust is the series filtered y
acf5=acf(out$y.robust,10, plot=FALSE)
plot(acf5, main="acf based on procedure (a)")
acf6=acf(out$y.robust,10,"partial", plot=FALSE)
plot(acf6, main="pcf based on procedure (a)")

# Procedure (B)
tank1=ARMAacf(ar = out$model$ar, lag.max = 10, pacf = FALSE)
tank2=ARMAacf(ar = out$model$ar, lag.max = 10, pacf = TRUE)
acf7=acf5
acf8=acf6
acf7$acf[,,1]=as.matrix(tank1)
plot(acf7, main="acf based on procedure (b)")
acf8$acf[,,1]=as.matrix(tank2)
plot(acf8, main="pcf based on procedure (b)")





