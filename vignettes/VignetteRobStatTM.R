## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
	keep.source=TRUE,
	tidy=TRUE,
	message=FALSE, 
	fig.path='Plots/', 
	fig.align='center', 
	fig.width=4, 
	fig.height=3, 
	fig.keep='last',
	fig.show='hide',
	dev.args=list(pointsize=10),
	tidy.opts=list(width.cutoff=40),
	cache=FALSE)
options(width=70)


## ----eval = FALSE,echo = T-----------------------------------------------
## install.packages("RobStatTM")


## ----echo = T------------------------------------------------------------
library("RobStatTM")


## ----eval = FALSE, echo = T----------------------------------------------
## system.file("scripts", package = "RobStatTM")


## ----echo = T------------------------------------------------------------
data(shock)
head(shock,2)


## ----echo = T------------------------------------------------------------
data(wood, package = "robustbase")
head(wood,1) 


## ----echo = T------------------------------------------------------------
minerall1 <- quantreg::rq(zinc ~ copper, data=mineral)


## ----echo = T------------------------------------------------------------
LSfit <- lm(zinc ~ copper, data = mineral)
control <- lmrobdet.control(family = "mopt",eff = 0.95)
robfit <- lmrobdetMM(zinc ~ copper, control = control, data = mineral)
fmLSrob <- fit.models(LSfit,robfit)
class(fmLSrob)
summary(fmLSrob)


## ----echo = T,eval = F,results = F---------------------------------------
## help(summary.lmfm)


## ----echo = T,eval = F,results = F---------------------------------------
## help(plot.lmfm)


## ----echo = T------------------------------------------------------------
args(plot.lmfm)


## ----echo = 2,results = F------------------------------------------------
png(file = "Plots/.png", width = 6, height = 4, units = "in",
   pointsize = 6, res = 600)
plot(fmLSrob,which.plots = 11)
dev.off()


## ----echo = 2,results = F------------------------------------------------
png(file = "Plots/qqplotResiduals.png", width = 6, height = 4, units = "in",
   pointsize = 6, res = 600)
plot(fmLSrob,which.plots = 2)
dev.off()


## ----echo = T,warning = F------------------------------------------------
library(robust)  # This is only needed until the package fit.models is updated in CRAN


## ----echo = T------------------------------------------------------------
data(wine)
wine5 <- wine[,1:5]
cov.fm <- fit.models(Classic = covClassic(wine5),Robust = covRob(wine5,type = "auto"))
class(cov.fm)
summary(cov.fm)


## ----echo = T------------------------------------------------------------
names(cov.fm)
names(cov.fm$Classic)
names(cov.fm$Robust)


## ----echo = T, eval = F--------------------------------------------------
## help(plot.covfm)


## ----echo = 2,results = F------------------------------------------------
png(file = "Plots/eigenvalues.png", width = 6, height = 4, units = "in",
   pointsize = 6, res = 600)
plot(cov.fm,which.plot = 2)
dev.off()


## ----echo = 2,results = F------------------------------------------------
png(file = "Plots/distances.png", width = 6, height = 4, units = "in",
   pointsize = 6, res = 600)
plot(cov.fm,which.plot = 3)
dev.off()


## ----echo = 2,results = F------------------------------------------------
png(file = "Plots/ellipses.png", width = 6, height = 4, units = "in",
   pointsize = 6, res = 600)
plot(cov.fm,which.plot = 4)
dev.off()


## ----echo = 2,results = F------------------------------------------------
png(file = "Plots/distanceDistance.png", width = 6, height = 4, units = "in",
   pointsize = 6, res = 600)
plot(cov.fm,which.plot = 5)
dev.off()

