# algae.R
# EXAMPLE 5.4
# Figures 5.14-5.15
# Table 5.4

library(RobStatTM)
data(algae)

#Robust fit
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

# We now recommend to use family "mopt" with efficiency = 0.95 as defaults.
# Using that as in next line results in almost no change in Figure 5.15
# cont <- lmrobdet.control(bb = 0.5, efficiency = 0.95, family = "mopt")
algaerob <- lmrobdetMM(V12 ~ ., data=algae, control=cont)

#LS fit
algaels <- lm(V12 ~ ., data=algae)

#LS fit without outliers
algaelsd <- lm(V12 ~ ., data=algae, subset= -c(36, 77))

#-----------------------------------------------------
#Fig 5.14
plot(algaels, which=2, pch=19)
abline(h=c(-2.5, 0, 2.5), lty=2)

#-------------------------------------------------------
#Fig 5.15
plot(algaerob, which=2, id.n=2, pch=19)
abline(h=c(-2.5, 0, 2.5)*algaerob$scale, lty=2)




