# mineral.R
# EXAMPLE 5.1  
# Figures 5.1 - 5.7
# Table 5.1

library(RobStatTM)

data(mineral)
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

#LS fit
mineralls <- lm(zinc ~  copper, data=mineral)

#L1 fit
minerall1 <- quantreg::rq(zinc ~  copper, data=mineral)

#LS without outlier fit
minerallssin <- lm(zinc ~ copper, data = mineral[-15,  ])

#Fig 5.1.
plot(zinc ~ copper, data=mineral, pch=19, cex=1.3)
abline(mineralls, lwd=2, col='red')
abline(minerall1, lwd=2, col='blue')
abline(minerallssin, lwd=2, col='green4')
text(c(600,600,600), c(29,55,82), c("LS-15","L1", "LS"), cex=1.3)
text(mineral[15,1], mineral[15,2]-6, "15", cex=1.3)



#Figure 5.2
plot(mineralls, which=2, add.smooth=FALSE, pch=19, id.n=2, cex.id = 1.2)
abline(h=c(2.5, 0, -2.5), lty=2, lwd=2)

#Figure 5.3
sigmaLS <- summary( mineralls)$sigma
plot(mineralls, which=1, add.smooth=FALSE, pch=19, id.n=2, cex.id = 1.2)
abline(h=c(2.5, 0, -2.5)*sigmaLS, lty=2, lwd=2)

#--------------------------------------------


#MM fit
mineralMM <- lmrobdetMM(zinc ~ copper, data = mineral, control=cont)

#Fig 5.4
plot(zinc ~ copper, data=mineral, pch=19, cex=1.3)
abline(minerallssin, lwd=2, col='green4')
abline(mineralMM, lwd=2, col='magenta')
text(c(600,600), c(45,29), c("ROB","LS-15"), cex=1.2)
text(mineral[15,1], mineral[15,2]-6, "15", cex=1.2)
#------------------------------------------------------------


#------------------------------------------------------------
#Fig 5.5
plot(mineralMM, which=4, add.smooth=FALSE, pch=19, cex.id=1.3)

#-----------------------------------------------

#Fig 5.6
plot(mineralMM, which=2, pch=19, id.n=3, cex.id=1.2)
abline(h=c(-2.5, 0, 2.5) * mineralMM$scale, lty=2, lwd=2)

#----------------------------------------------------------

#Fig 5.7
plot(sort(abs(resid(mineralls)))[-53], sort(abs(resid(mineralMM)))[-53], xlab="Least Squares residuals",
     ylab="Robust residuals", pch=19, cex=1.3)
abline(0, 1, lwd=2, col='red')
















