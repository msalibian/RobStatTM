 
library(RobStatTM)
#data
u=1:5
block=rep(u,8)
block=as.factor(block)
variety=rep(1,5)
for ( i in 2:8) 
variety=c(variety,rep(i,5))
variety=as.factor(variety)
response=scan("oats.txt")
# response = data(oats)
response1=response[1:40]
response2=response[41:80]
 
  
##LS regression
 oats1LS=lm(response1~variety+block )
oats1LS_var=lm(response1~ block )
 oats1LS_block=lm(response1~variety )
oats2LS=lm(response2~variety+block )
oats2LS_var=lm(response2~ block )
oats2LS_block=lm(response2~variety  )
sigma1LS=summary(oats1LS)$sigma
sigma2LS=summary(oats2LS)$sigma
cont=lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")





##Classical ANOVA  tests
anov1_var =anova(oats1LS,oats1LS_var)
anov1_block=anova(oats1LS,oats1LS_block)
anov2_var =anova(oats2LS,oats2LS_var)
anov2_block=anova(oats2LS,oats2LS_block)


 
##M regressions
oats1M=lmrobM(response1~variety+block,control=cont )
oats1M_var=lmrobM(response1~ block, control=cont )
 oats1M_block=lmrobM(response1~variety,control=cont )
oats2M=lmrobM(response2~variety+block,control=cont )
oats2M_var=lmrobM(response2~ block,control=cont )
oats2M_block=lmrobM(response2~variety,control=cont )
sM2=oats2M$scale
 
##Robust ANOVA  tests

anov1M_var =rob.linear.test(oats1M,oats1M_var)
anov1M_block=rob.linear.test(oats1M,oats1M_block)
anov2M_var =rob.linear.test(oats2M,oats2M_var)
anov2M_block=rob.linear.test(oats2M,oats2M_block)
 
# qqplot of LS residuals of modified data   (Fig 4.2)
#pdf(file='Figure 4-3.pdf', bg='transparent', width=10, height=8)
qqnorm(oats2LS$resid/sigma2LS,ylab="residuals")
lines(c(-2.5,2.5),c(0,0),lty=2)
#dev.off()

 

# qqplot of robust  residuals of modified data  (Fig 4.4)
u=oats2M$resid[c(24,36,1,35,20)]/sM2 
w=c("36","24","1","35","20")
v=c(qnorm(.5/40)+.05,qnorm(1.5/40)+0.05,qnorm(37.5/40)+.05,qnorm(38.5/40)+.05,qnorm(39.5/40)+.05)
#pdf(file='Figure 4-4.pdf', bg='transparent', width=10, height=8)
qqnorm(oats2M$resid/sM2,ylab="residuals")
lines(c(-3,3),c(2.5,2.5),lty=2)
lines(c(-3,3),c(-2.5,-2.5),lty=2)
lines(c(-3,3),c(0,0),lty=2)
text(v,u,w)
#dev.off()

 

 
 
 



