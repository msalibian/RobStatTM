# wood.R
# EXAMPLE 5.2
# Figures 5.8 - 5.12
# wood data from robustbase

# N.B. The 4 figures produced by the code below are very similar to Figures
# 5.8 - 5.12 in the book, except the vertical scales are different due to
# the use of the plot method

library(RobStatTM)
data(wood, package='robustbase')
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

# We now recommend to use family "mopt" with efficiency = .95 as defaults
# Using those  results in almost no change in Figures 5.11-5.12

#MM fit
woodMM <- lmrobdetMM(y ~ ., data=wood, control=cont)

#LS fit
woodLS <- lm(y ~ ., data=wood)

#-------------------------------------------------------
#Fig 5.8
# Nothing happens in this figure
# text(woodLS$fitted[28]+.03,woodLSresst[28],"28")
sigmaLS <- summary(woodLS)$sigma
plot(woodLS, which=1, add.smooth=FALSE, pch=19, id.n=2, cex.id = 1.2)
abline(h=c(-2.5, 0, 2.5) * sigmaLS, lty=2, lwd=2)

#----------------------------
#Figure 5.9
plot(woodLS, which=2, pch=19)


#-------------------------------------------
# Fig. 5.10
plot(woodMM, which=4, add.smooth=FALSE, pch=19, cex.id=1.3, id.n=4)


#--------------------------
#Figure 5.11
plot(woodMM, which=2, add.smooth=FALSE, pch=19, cex.id=1.3, id.n=4)
abline(h=c(-2.5, 0, 2.5) * woodMM$scale, lty=2)
#----------------------------

#Figure 5.12
wq <- which( abs(resid(woodMM)) > 2.5 * woodMM$scale )
plot(sort(abs(resid(woodLS)[-wq] ) ), sort(abs(resid(woodMM)[-wq])), pch=19, xlim=c(0, .05), ylim=c(0, .05))
abline(0,1)












