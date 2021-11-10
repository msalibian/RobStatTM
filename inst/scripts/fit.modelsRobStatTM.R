# Section 1: Introduction
install.packages("fit.models")
library(RobStatTM)
library("fit.models")

# Section 2: fit.models for LS and lmrobdetMM model fits
fmclass.add.class("lmfm","lmrobdetMM")
LSfit <- lm(zinc ~ copper, data = mineral)
control <- lmrobdet.control(family = "mopt",eff = 0.95)
# The choices "mopt" and 0.95 are defaults
robfit <- lmrobdetMM(zinc ~ copper, control = control, data = mineral)
fmLSrob <- fit.models(LSfit,robfit)
class(fmLSrob)
names(fmLSrob)
round(coef(fmLSrob),3)
summary(fmLSrob)

help(summary.lmfm)
help(plot.lmfm)

plot(fmLSrob) # View plot types sequentially at the RStudio Console
plot(fmLSrob,which.plots = "ask") # Choose plot types at Console
plot(fmLSrob,which.plots = 10)
plot(fmLSrob,which.plots = 7)
plot(fmLSrob,which.plots = 1)
plot(fmLSrob, which.plots = c(3,4)) # View two plot types at Console

fmLSonly <- fit.models(LSfit)
class(fmLSonly)
names(fmLSonly)
coef(fmLSonly)
plot(fmLSonly, which.plots = "ask")

# Section 3: fit.models for covClassic and covRob
data(wine)
wine3 <- wine[,1:3]
fmCovRob <- fit.models(Classic = covClassic(wine3),
                       Robust = covRob(wine3,type = "auto"))
class(fmCovRob)
summary(fmCovRob)

names(fmCovRob)
names(fmCovRob$Classic)
names(fmCovRob$Robust)

help(plot.covfm)

plot(fmCovRob,which.plot = 1)
plot(fmCovRob,which.plot = 2)
plot(fmCovRob,which.plot = 3)
plot(fmCovRob,which.plot = 4)


