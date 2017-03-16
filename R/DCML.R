
# require(lmrob)
# require(penseinit)
# require(quantreg)
 

mscale <- function(u, tol, delta=.5, max.it=100, tuning.chi=1.5477) {
  # M-scale of a sample u  
  # tol: accuracy
  # delta: breakdown point (right side)
  # Initial 
  s0 <- median(abs(u))/.6745
  err <- tol + 1
  it <- 0
  while( (err > tol) & ( it < max.it) ) {
    it <- it+1
    s1 <- sqrt( s0^2 * mean(rho((u/s0), tuning.chi)) / delta )
    err <- abs(s1-s0)/s0
    s0 <- s1
  }
  return(s0)
}



covdcml=function(resLS,res0,C,sig0,t0,p,n)
##Computation of the asymptotic covariance matrix of the DCML estimator
{t0=1-t0
c=mean(psibs(res0/sig0,3.44)*resLS)
a1=mean(psibs(res0/sig0,3.44)^2)
b=mean(psibspri(res0/sig0,3.44))
deltasca=0.5*(1-(p/n))
a2=mscale(resLS,.00001,deltasca)^2
tuti=t0^2*sig0^2*a1/b^2 + a2*(1-t0)^2 +2*t0*(1-t0)*sig0*c/b
V=tuti*solve(C)
V
}


rho <- function(u, cc=1.5477) {
  w <- as.numeric( abs(u) <= cc )
  v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
  v <- v*6/cc^2
  return(v)
}

rhoint <- function(e)
  return(integrate(function(a, cc) rho(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value)


find.tuning.chi <- function(delta, low=.5, upp=10) {
  return( uniroot( function(e) (rhoint(e)-delta), lower=low, upper=upp)$root )
}


rhoprime <- function(r, cc) { 
  # bisquare rho' = psi function
  r <- r / cc
  gg <- r*(1-r^2)^2*as.numeric(abs(r)<=1)
  return(gg)
}

rhoprime2 <- function(r, cc) {
  #Derivative of the bisquare psi function
  r <- r/cc
  gg <- (1-r^2)*(1-5*r^2)*as.numeric(abs(r)<=1)/cc
  return(gg)
}

# effi <- function(e) {
#   a <- integrate(function(a, cc) (rhoprime(a, cc)^2)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
#   b <- integrate(function(a, cc) rhoprime2(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
#   return( 1/(a/b^2) ) # efficiency!
# }
# 
# uniroot( function(e) (effi(e)-.85), lower=3, upper=4)$root

  psibs=function(t,c)  
# bisquare psi function
 {r=t/c
 gg=r*(1-r^2)^2*(abs(r)<=1)
 gg
}
        
 
 psibspri=function (t,c)   
#Derivative of the bisquare psi function
{ r=t/c
 gg=(1-r^2)*(1-5*r^2)*(abs(r)<=1)/c
gg
}
        
 MMPY <- function(X, y, control, mf, corr.b=control$corr.b) {
   # This function will be called from lmrob, so control will be valid
   # X will already contain a column of ones if needed
   #Compute an MM-estimator taking as initial Pe?a Yohai
   #INPUT
   #X nxp matrix, where n is the number of observations and p the number of  columns
   #y vector of dimension  n with the responses
   #
   #OUTPUT
   #outMM output of the MM estimator (lmrob) with 85% of efficiency and PY as initial
   #
   # 
   # cont1=lmrob.control(tuning.psi=3.44)
   n <- nrow(X)
   p <- ncol(X)
   # if(intercept==TRUE )
   # { p=p+1}
   dee <- control$bb
   if(corr.b) dee <- dee * (1-(p/n))
   a <- pyinit(X=X, y=y, intercept=FALSE, deltaesc=dee, 
               cc.scale=control$tuning.chi, 
               prosac=control$prosac, clean.method=control$clean.method, 
               C.res = control$C.res, prop=control$prop, 
               py.nit = control$py.nit, en.tol=control$en.tol, 
               mscale.maxit = control$mscale.maxit, mscale.tol = control$mscale.tol,
               mscale.rho.fun=control$mscale.rho.fun)
   # betapy2 <- a$initCoef[,1]
   # sspy2 <- a$objF[1]
   # refine the PY candidates to get something closer to an S-estimator for y ~ X1
   kk <- dim(a$initCoef)[2]
   best.ss <- +Inf
   for(i in 1:kk) {
     tmp <- refine.sm(x=X, y=y, initial.beta=a$initCoef[,i], 
                      #initial.scale=initial$objF[1], 
                      k=control$refine.PY, conv=1, b=dee, cc=control$tuning.chi, step='S')
     if(tmp$scale.rw < best.ss) {
       best.ss <- tmp$scale.rw # initial$objF[1]
       betapy <- tmp$beta.rw # initial$initCoef[,1]
     }
   }
   # uu <- list(coef=betapy, scale=best.ss)
   # betapy2 <- betapy
   # sspy2 <- best.ss
   S.init <- list(coef=betapy, scale=best.ss)
   control$method <- 'M'
   control$cov <- ".vcov.w"
   # # lmrob() does the above when is.list(init)==TRUE, in particular:
   #
   # if (control$method == "MM" || substr(control$method, 1, 1) == "S")
   #   control$method <- substring(control$method, 2)
   # ## check for control$cov argument
   # if (class(init)[1] != "lmrob.S" && control$cov == '.vcov.avar1')
   #   control$cov <- ".vcov.w"
   
   outMM <- lmrob.fit(X, y, control, init=S.init, mf=mf)
   # if(intercept==TRUE)
   #    { outMM=lmrob(y~X, control=cont1,init=uu)}else
   #    { outMM=lmrob(y~X-1, control=cont1,init=uu)}
   return(outMM)
 }

old.MMPY <- function(X, y, intercept=FALSE) {
   # X will already contain a column of ones if needed
   #Compute an MM-estimator taking as initial Pe?a Yohai
   #INPUT
   #X nxp matrix, where n is the number of observations and p the number of  columns
   #y vector of dimension  n with the responses
   #
   #OUTPUT
   #outMM output of the MM estimator (lmrob) with 85% of efficiency and PY as initial
   #
   # 
  
   cont1 <- lmrob.control(tuning.chi = 1.5477, bb = 0.5, tuning.psi = 3.4434)
   n <- nrow(X)
   p <- ncol(X)
   if(intercept)
   { p=p+1}
   dee <- .5*(1-(p/n))
   a <- pyinit(X=X, y=y, intercept=FALSE, deltaesc=dee, 
               cc.scale=cont1$tuning.chi, 
               prosac=.5, clean.method='threshold', C.res = 2, prop=.2, 
               py.nit = 20, en.tol=1e-5, mscale.rho.fun='bisquare')
   betapy2 <- a$initCoef[,1]
   sspy2 <- a$objF[1]
   S.init <- list(coef=betapy2, scale=sspy2)
   # print('calling lmrob.fit')
   # control$method <- 'M'
   # outMM <- lmrob.fit(X, y, control, init=S.init, mf=mf)
   if(intercept)
      { outMM=lmrob(y~X, control=cont1,init=S.init)}else
      { outMM=lmrob(y~X-1, control=cont1,init=S.init)}
   return(outMM)
 }
 
  

DCML_FINAL=function(X,y, outMM, intercept=TRUE)
#INPUT
#X nxp matrix, where n is the number of observations and p the number of  columns, the 1's of the intercept are not included
#y vector of dimension  n with the responses
#intercept FALSE the regression does not includes intercept, TRUE, the regression includes intercept
#outMM output of MMPY or SM_PY
#OUTPUT
#coef vector of coefficients, the first element is the intercept when there is one
#cov covariance matrix of the coefficients
#resid vector of residuals
#weight  vector with weights that the MMestimator assigns to every observation
#sigma standad error of errot term
 

{
XX=X
n=nrow(X)
p=ncol(X)
if(intercept==TRUE )
    {XX=cbind(rep(1,n),X)
    p=p+1}
cont1=lmrob.control(tuning.psi=3.44)
if(intercept==TRUE)
   {out3=lm(y~X)}else
   {out3=lm(y~X-1)}
betaLS=out3$coeff
resLS=out3$resid
dee=.5*(1-(p/n))
beta0=outMM$coeff
weight=outMM$rweights
residuos=outMM$resid
sigma=mscale(residuos,.00001,dee)
##Begin the computation of the DCML
deltas=.3*p/n 
C=t(XX)%*%diag(weight)%*%XX/sum(weight) 	
d=(beta0-betaLS)%*%C%*%(beta0-betaLS)
d=d/sigma^2
t0=min(1,sqrt(deltas/d))
beta1=t0*betaLS+(1-t0)*beta0
V=covdcml (resLS,residuos,C,sigma,t0,p,n)/n	
resid=y-XX%*%beta1
sigma=mscale(resid, .00001,dee) 
out=list(coef=beta1, cov=V, resid=resid,   sigma=sigma )
out
}

SMPY <- function(mf, y, control, split, corr.b=control$corr.b) { 
  if(missing(control)) 
    control <- lmrob2.control(tuning.chi = 1.5477, bb = 0.5, tuning.psi = 3.4434)
  int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1)
  if(missing(split)) {
    split <- splitFrame(mf, type=control$split.type) 
  }
  Z <- split$x1 # x1 = factors, x2 = continuous, if there's an intercept it's in x1!
  X <- split$x2
  if(int.present) Z <- Z[, -1]
  n <- nrow(X)
  q <- ncol(Z)
  p <- ncol(X)
  dee <- control$bb
  gamma <- matrix(NA, q, p)
  #Eliminate Z from X by L1 regression , obtaining a matrix X1
  oldw <- options()$warn
  options(warn=-1)
  for( i in 1:p) gamma[,i] <- coef(rq(X[,i]~Z-1))
  options(warn=oldw)
  # Now regress y on X1, find PY candidates
  X1 <- X - Z %*% gamma
  pp <- p + if(int.present) 1 else 0
  if(corr.b) dee <- dee*(1-(pp/n))
  initial <- pyinit(intercept=int.present, X=X1, y=y, 
                    deltaesc=dee, cc.scale=control$tuning.chi, 
                    prosac=control$prosac, clean.method=control$clean.method, 
                    C.res = control$C.res, prop=control$prop, 
                    py.nit = control$py.nit, en.tol=control$en.tol, 
                    mscale.maxit = control$mscale.maxit, mscale.tol = control$mscale.tol,
                    mscale.rho.fun=control$mscale.rho.fun)
  # beta <- uu$coef  # refine the PY candidates to get something closer to an S-estimator for y ~ X1
  kk <- dim(initial$initCoef)[2]
  best.ss <- +Inf
  Xtmp <- X1
  if(int.present) Xtmp <- cbind(rep(1, nrow(Xtmp)), Xtmp)
  for(i in 1:kk) {
    tmp <- refine.sm(x=Xtmp, y=y, initial.beta=initial$initCoef[,i], 
                     #initial.scale=initial$objF[1], 
                     k=control$refine.PY, conv=1, b=dee, cc=control$tuning.chi, step='S')
    if(tmp$scale.rw < best.ss) {
      best.ss <- sspy <- tmp$scale.rw # initial$objF[1]
      betapy <- tmp$beta.rw # initial$initCoef[,1]
    }
  }
  uu <- list(coef=betapy, scale=sspy)

  # y1 <- as.vector(y - X1 %*% beta)
  #
  # Do an M-step starting from this S?
  #
  if(int.present) {
    out0 <- lmrob(y~X1, control=control,init=uu) } else {
      out0 <- lmrob(y~X1 - 1, control=control,init=uu)
    }
  beta <- out0$coeff
  y1 <- resid( out0 )
  #
  # after eliminating the influence of X1 we make an L1 regression using as covariables Z
  # fit the model I( y1 - X1 %*% hat(beta) ) ~ Z
  options(warn=-1)
  fi <- coef(rq(y1~Z-1))
  options(warn=oldw)
  # retransform the coefficients to the original matrices X and Z
  if(int.present) {
    beta00 <- c(beta, fi - gamma %*% beta[-1])
  } else { beta00 <- c(beta, fi - gamma %*% beta) }
  res <- as.vector( y1 - Z%*%fi )
  # XX <- cbind(X,Z)
  nc <- ncol(X) + ncol(Z) + if(int.present) 1 else 0
  dee <- control$bb
  if(corr.b) dee <- dee * (1-(nc/n))
  ss <- mscale(u=res, tol=1e-7, delta=dee, tuning.chi=control$tuning.chi)
  XX <- model.matrix(attr(mf, 'terms'), mf)
  uu <- list(coef=beta00, scale=ss, residuals=res)
  # print(uu)
  # tmp <- refine.sm(x=XX, y=y, initial.beta=beta00,
  #                   k=500, conv=1, b=dee, cc=control$tuning.chi, step='S')
  # print(tmp$conver)
  # uu <- list(coef=as.vector(tmp$beta.rw), scale=tmp$scale.rw)
  # print(uu)
  # Compute the MMestimator using lmrob
  control$method <- 'M' 
  control$cov <- ".vcov.w"
  # lmrob() sets the above when is.list(init)==TRUE
  #
  outlmrob <- lmrob.fit(XX, y, control, init=uu, mf=mf)
  return(c(outlmrob, init.SMPY=uu)) 
}

old.SMPY=function(y,X,Z, intercept=TRUE)
   #old penayohai2 Compute an MM-estimator for mixed model using lmrob and taking as initial an SM estimator based on Pe?a-Yohai
   #INPUT
   #y response vector
   #X matrix of continuous covariables
   #Z matrix of  qualitative covariables
   #OUTPUT
   #out_lmrob output of lmrob take as initial a S-M estimator for mixed models computed with Pe?a Yohai
 {
    cont1=lmrob.control(tuning.psi=3.44) 
    n=nrow(X)
    q=ncol(Z)
    p=ncol(X)
    hh1=1+intercept
    hh2=p+intercept
    hh3=q+intercept
    gamma=matrix(0, (q+intercept),p)
    
    
    for( i in 1:p)
    {gamma[,i]=rq(X[,i]~Z)$coeff}
    X1=X-Z%*%gamma[hh1:hh3,]
    #aa=c(dim(Z),dim(X),dim(X1),dim(y))
    
    
    dee=.5*(1-((p+1)/n))
    out= pyinit(intercept=intercept,X=X1, y=y, deltaesc=dee, cc.scale=1.547, prosac=.5,clean.method="threshold", C.res = 2, prop=.2, py.nit = 20, en.tol=1e-5)
    betapy=out$initCoef[,1]
    sspy=out$objF[1]
    uu=list(coeff=betapy,scale=sspy)
    out0=lmrob(y~X1, control=cont1,init=uu) 
    beta=out0$coeff
    
    y1=y-X1%*%beta[hh1:hh2]
    
    
    fi=rq(y1~Z)$coeff
    oo=NULL
    if(intercept==TRUE)
      oo=fi[1]
    #print(length(beta))
    #print(dim(gamma[hh1:hh3,]))
    tt=gamma[hh1:hh3,]
    if(p==1)
      tt=matrix(tt,q,p)
    beta00=c(oo,beta[hh1:hh2],fi[hh1:hh3]-  tt%*%beta[hh1:hh2])
    res=y1-fi[1]-Z%*%fi[hh1:hh3]
    res=as.vector(res)
    dee=.5*(1-((p+q+intercept)/n))
    ss=mscale(res,.0001,dee)
    uu=list(coeff=beta00,scale=ss)
    XX=cbind(X,Z)
    beta0=lmrob(y~XX,control=cont1,init=uu)$coeff
    beta0
 }
 

# splitFrame <- function(mf, x = model.matrix(mt, mf),
#                        type = c("f","fi", "fii"))
# {
#   mt <- attr(mf, "terms")
#   type <- match.arg(type)
#   x <- as.matrix(x)
#   p <- ncol(x)
#   
#   ## --- split categorical and interactions of categorical vars.
#   ##     from continuous variables
#   factors <- attr(mt, "factors")
#   factor.idx <- attr(mt, "dataClasses") == "factor"
#   if (!any(factor.idx)) ## There are no factors
#     return(list(x1.idx = rep.int(FALSE, p), x1=matrix(,nrow(x),0L), x2=x))
#   switch(type,
#          ## --- include interactions cat * cont in x1:
#          fi = { factor.asgn <- which(factor.idx %*% factors > 0) },
#          ## --- include also continuous variables that interact with factors in x1:
#          ##     make sure to include interactions of continuous variables
#          ##     interacting with categorical variables, too
#          fii = { factor.asgn <- numeric(0)
#          factors.cat <- factors
#          factors.cat[factors.cat > 0] <- 1L ## fix triple+ interactions
#          factors.cat[, factor.idx %*% factors == 0] <- 0L
#          for (i in 1:ncol(factors)) {
#            comp <- factors[,i] > 0
#            ## if any of the components is a factor: include in x1 and continue
#            if (any(factor.idx[comp])) {
#              factor.asgn <- c(factor.asgn, i)
#            } else {
#              ## if there is an interaction of this term with a categorical var.
#              tmp <- colSums(factors.cat[comp,,drop=FALSE]) >= sum(comp)
#              if (any(tmp)) {
#                ## if no other continuous variables are involved
#                ## include in x1 and continue
#                ## if (identical(factors[!comp, tmp], factors.cat[!comp, tmp]))
#                if (!all(colSums(factors[!factor.idx & !comp, tmp, drop=FALSE]) > 0))
#                  factor.asgn <- c(factor.asgn, i)
#              }
#            }
#          } },
#          ## --- do not include interactions cat * cont in x1:
#          f = { factor.asgn <- which(factor.idx %*% factors & !(!factor.idx) %*% factors) },
#          stop("unknown split type"))
#   x1.idx <- attr(x, "assign") %in% c(0, factor.asgn) ## also include intercept
#   names(x1.idx) <- colnames(x)
#   
#   ## x1: factors and (depending on type) interactions of / with factors
#   ## x2: continuous variables
#   list(x1 = x[,  x1.idx, drop=FALSE],
#        x1.idx = x1.idx,
#        x2 = x[, !x1.idx, drop=FALSE])
# }

