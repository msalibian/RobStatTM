# ar3.R
# Simulated AR3 data
# true pararameters 4/3, -5/6 , 1/6
# Table 8.1
 
library(RobStatTM)
# Must install robustarima
library(robustarima)


set.seed(600)
n.innov = 300
n = 200
phi=c(4/3, -5/6,1/6  )

n.start = n.innov - n 
innov = rnorm(n.innov)
cont=lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")
ar3= arima.sim(model = list(ar = phi), n, innov = innov, n.start = n.start)

 

 
#no outliers
ar3lm=lm(ar3[4:200]~ ar3[3:199]+ar3[2:198]+ar3[1:197])
ar3mm=lmrobdetMM(ar3[4:200]~ ar3[3:199]+ar3[2:198]+ar3[1:197] )
ar3tau=arima.rob(ar3~1, p=3) 
ar3tau$regcoef=ar3tau$regcoef*(1-sum(ar3tau$model$ar))
 
summary(ar3lm)
summary(ar3mm)
ar3tau

#5% of outliers
 
  
ao=rep(0,n)
tt=seq(20,200,20)
ao[tt]=4
 
 

ar3c5 =ar3  +ao
ar3c5lm=lm(ar3c5[4:200]~ ar3c5[3:199]+ar3c5[2:198]+ar3c5[1:197])
ar3c5mm=lmrobdetMM(ar3c5[4:200]~ ar3c5[3:199]+ar3c5[2:198]+ar3c5[1:197], control=cont )
ar3c5tau=arima.rob(ar3c5~1, p=3) 
ar3c5tau$regcoef=ar3c5tau$regcoef*(1-sum(ar3c5tau$model$ar))
summary(ar3c5lm)
summary(ar3c5mm)
ar3c5tau

 
#10% outliers
 

 ao=rep(0,n)
tt=seq(10,200,10)
ao[tt]=4
	
ar3c10=ar3 +ao
ar3c10lm=lm(ar3c10[4:200]~ ar3c10[3:199]+ar3c10[2:198]+ar3c10[1:197])
ar3c10mm=lmrobdetMM(ar3c10[4:200]~ ar3c10[3:199]+ar3c10[2:198]+ar3c10[1:197],control=cont)
ar3c10tau=arima.rob(ar3c10~1, p=3) 
ar3c10tau$regcoef=ar3c5tau$regcoef*(1-sum(ar3c5tau$model$ar))
summary(ar3c10lm)
summary(ar3c10mm)
ar3c10tau





 
 