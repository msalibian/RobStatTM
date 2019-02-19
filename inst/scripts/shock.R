# shock.R
# EXAMPLE 4.1
# Figures 4.1
# Table 4.1

library(RobStatTM)

data(shock)
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

#LS fit
shockls <- lm(time ~ n.shocks, data = shock)

#LS fit without outliers
shockls124 <- lm(time ~ n.shocks, data = shock, subset = -c(1, 2, 4))

#--------------------
#Figure 4.1
plot(time ~ n.shocks, data=shock, xlab="number of shocks", ylab="average time", pch=19, cex=1.3)
abline(shockls124, lwd=2, col='gray30')
abline(shockls, lwd=2, col='tomato')
text(shock[c(1,2,4),1],shock[c(1,2,4),2]-.4, labels=c("1","2","4"), cex=1.2)
text(1,10.5,"LS", col='tomato', cex=1.3)
text(1,7.2,"LS -", col='gray30', cex=1.3)
#---------------------------------

#L1 fit
shockl1 <- quantreg::rq(time~n.shocks , data=shock)

# M fit
shockrob <- lmrobM(time ~ n.shocks, data = shock,control=cont)


#--------------------------
#Figure 4.3
plot(time ~ n.shocks, data=shock, xlab="number of shocks", ylab="average time", pch=19, cex=1.3)
abline(shockls124, lwd=2, col='gray30')
abline(shockls, lwd=2, col='tomato')
abline(shockrob, lwd=2, col='green')
abline(shockl1, lwd=2, col='blue')
text(shock[c(1,2,4),1],shock[c(1,2,4),2]-.4, labels=c("1","2","4"), cex=1.3)
text(1,10.5,"LS", cex=1.3)
text(1,7.2,"LS -", cex=1.3)
text(1,8.3,"L1", cex=1.3)
text(1,7.7,"M", cex=1.3)














