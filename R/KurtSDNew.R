#' Robust multivariate location and scatter estimators
#'
#' This function computes robust multivariate location and scatter
#' estimators using both random and deterministic starting points.
#'
#' This function computes robust multivariate location and scatter
#' using both Pen~a-Prieto and random candidates.
#'
#' @export KurtSDNew initPP
#' @aliases KurtSDNew initPP
#' @rdname KurtSDNew
#'
#' @param X a data matrix with observations in rows.
#' @param muldirand used to determine the number of random directions (candidates), which
#' is \code{max(p*muldirand, dirmin)}, where \code{p} is the number of columns in \code{X}.
#' @param muldifix used to determine the number of random directions (candidates), which
#' is \code{min(n, 2*muldifix*p)}.
#' @param dirmin minimum number of random directions
#'
#' @return A list with the following components:
#' \item{idx}{A zero/one vector with ones in the positions of the suspected outliers}
#' \item{disma}{Robust squared Mahalanobis distances}
#' \item{center}{Robust mean estimate}
#' \item{cova}{Robust covariance matrix estimate}
#' \item{t}{Outlyingness of data points}
#'
#' @author Ricardo Maronna, \email{rmaronna@retina.ar}, based on original code
#' by D. Pen~a and J. Prieto
#'
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(bus)
#' X0 <- as.matrix(bus)
#' X1 <- X0[,-9]
#' tmp <- initPP(X1)
#' round(tmp$cov[1:10, 1:10], 3)
#' tmp$center
#'
initPP <- KurtSDNew <- function(X, muldirand=20, muldifix=10,dirmin=1000) {

  oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()))

n0=dim(X)[1]; p=dim(X)[2]
nmxpass = 10         # max number of passes through basic outlier identification
#nmxpass = 1;         # max number of passes through basic outlier identification

# Cutoff values
ctf = 3.05 + 0.081*p
# minimum numbers of observations that are acceptable as good observations
lmn1 =max(c(floor(0.5*(n0 + p + 1)),2*p))
lmn2 = min(c(2*p,max(p + 2,floor(0.25*n0))))

#standardize data
mm = colMeans(X); S = cov(X)
xn =scale(X,center=mm,scale=FALSE)
Rr = chol(S)
xn= t(solve(t(Rr),t(xn)))
xs=xn
xn0 = cbind(xn, 1:n0)
x = xn0

# Initialize parameters
n = n0;
flg = 1
ndir = 0; pass = 1

while ( (n  > lmn1) && (pass <= nmxpass) && (flg == 1)) {
# Compute fixed directions
  mm=p*muldifix; XN=x[,1:p]
  sig=sqrt(colSums(t(XN)^2)); V=scale(t(XN),scale=sig, center=F)
  if (mm<floor(n/2) ) {
      ind=order(sig); V=V[,ind]
  #      		print(c(n0,n, mm, dim(V)))
      V=cbind(V[,1:mm], V[,(n-mm+1):n])
  }
  V1=V
  #Compute random directions
  set.seed(111)   #to get reproducible results
  nd =muldirand*p; nd=max(c(nd,dirmin))
  ndk = c(1,1,nd)
  V2 = KurNwm(XN,ndk)
  Vv=cbind(V1,V2)     #matrix of directions
  Pro = XN%*%Vv        #Projected data
  mupr = apply(Pro,2,median)
  prcen=scale(Pro,center=mupr,scale=FALSE)
  prcen=abs(prcen)
  madpr=apply(prcen,2,median)
  tt = scale(prcen,center=FALSE, scale=madpr)
  t =apply(tt/ctf,1,max)
  # Exit if no observations were labelled as outliers
  if (sum(t > 1)==0) {flg=0; break}

  # Exit if too many observations were labelled as outliers
  inn=(t<=1); nu = sum(inn)
  if (nu <= lmn2) {
        ix = order(t); ixx = ix[1:lmn2]
        x = x[ixx,]
        n=dim(x)[1]; pp=dim(x)[2]
        flg = 0;   break
  } #end if

## Extract remaining observations left to be analyzed and redo
  x = x[inn,];
  n=dim(x)[1]; pp=dim(x)[2]
  pass = pass + 1;
} #end while

#Extracting the indices of the outliers
idx = rep(1,n0)
nn=dim(x)[1]; pp=dim(x)[2]
idx[x[,pp]] = rep(0,nn);
    #recheck observations and relabel them
    #Cutoffs for Mahalanobis distances
    ctf = qchisq(.99,p)
#Mahalanobis distances using center and scale estimators based on good observations
  sidx = sum(idx)
  s1 = which(idx>0);   s2 = which(idx == 0)
 if (sidx > 0) {xx1 = xn0[s1,]; xx1r = xn[s1,]   }
 xx2r = xs[s2,]
 mm = colMeans(xx2r);   Ss = cov(xx2r)
if (sidx > 0) {dms=mahalanobis(xx1r,mm,Ss)
} else {dms=NULL}

#Ensure that at least lmn1 observations are considered good
ado = sidx + lmn1 - n0
if (ado > 0) { idms = order(dms); idm1 = idms[1:ado]
  s3 = s1[idm1];
  idx[s3] = rep(0,ado)
  sidx = sum(idx); s1 = (idx>0); s2 = (idx == 0);
  if (sidx > 0) {xx1 = xn0[s1,]; xx1r = xn[s1,] }
  xx2r = xs[s2,]
  mm = colMeans(xx2r); Ss = cov(xx2r)
  if (sidx > 0) {dms=mahalanobis(xx1r,mm,Ss) }
}
  # Check remaining observations and relabel if appropriate

 while (sidx > 0) {
  s1 =(dms <= ctf); s2 = sum(s1)
   #   print(c("s2", s2))
  if (s2 == 0) break
  # print(c("s1 pp",length(s1), pp  ))
      xx1=as.matrix(xx1);
    #  print(c("xx1, pp,s1",  dim(xx1), pp,length(s1)))
      if (length(s1)==1) xx1=t(xx1)
      if (pp>dim(xx1)[2]) print(c("xx1 ,pp", dim(xx1),pp))
      s3 =xx1[s1,pp]
  idx[s3] = rep(0,s2)
   sidx = sum(idx);
   s1 = (idx>0); s2 = (idx == 0);
  if (sidx > 0) {xx1 = xn0[s1,]; xx1r = xn[s1,];   }
   xx2r = xn[s2,]
  mm = colMeans(xx2r); Ss = cov(xx2r);
  if (sidx > 0) dms=mahalanobis(xx1r,mm,Ss)
}  #end while

# Values to be returned
s2 = (idx == 0); mm = colMeans(X[s2,]); Ss = cov(X[s2,]);
dms = mahalanobis(X,mm,Ss)
resu=list(idx=idx, disma=dms, center=mm, cova=Ss, t=t)
return(resu)
}

######## Auxiliary functions
KurNwm <- function(xs,ndk) {
   #
   # Computation of directions that maximize and minimize
   # the kurtosis coefficient of the projections
   # and other directions based on ordering the observations
   # The values of the projections are also computed
   # To be used as a subroutine of KurMain
   #
   # Inputs:  x, observations (by rows)
   #          ndir,  [ num_dir_max_kur num_dir_min_kur num_dir_sd ]
   # Outputs: V, directions maximizing/minimizing kurtosis coef.
   #             maximization directions first
   #

   # Initialization of parameters

   tol = 1.0e-10
   stdsel = 1

   sx0 = dim(xs)
   n = sx0[1]
   p0 = sx0[2]

   dkmax = 2      # Generate both minimization and maximization directions (do not change lightly)

   ## Initialization of vectors

   Vv = matrix(0,p0,0)

   ## Standardize data

   if (stdsel) {
      en = matrix(1,n,1)
      mm = apply(xs,2,mean)
      S = cov(xs)
      x0 = xs - en %*% mm
      Rr = chol(S)
      x0 = t(solve(t(Rr),t(x0)))
   } else x0 = xs

   # Computing directions

   ## Choice of minimization/maximization/others

   dk = 1

   ## Main loop to compute the desired number of directions

   while (dk <= dkmax) {
      xx = x0
      p = p0
      pmx = min(ndk[dk],p-1)
      M = diag(p0)

      for (i in 1:pmx) {
         if (dk == 1) av = MaxKur(xx) else if (dk == 2) av = MinKur(xx)

         a = av$v

         la = length(a)
         za = matrix(0,la,1)
         za[1] = 1
         w = a - za
         nw = t(w) %*% a
         if (abs(nw) > tol) Q = diag(la) - (w %*% t(w))/c(nw) else Q = diag(la)

         ## Compute projected values

         if ((i > 1)||(dk > 1)) Vv = cbind(Vv,M %*% a) else Vv = M %*% a

         Qp = Q[,2:p]
         M = M %*% Qp

         ## Reduce dimension

         Tt = xx %*% Q
         xx = Tt[,2:p]

         p = p - 1
      }

      ## Compute last projection (if needed)

      if (ndk[dk] > (p0 - 1)) Vv = cbind(Vv,M)

      ## Proceed to next class of directions

      dk = dk + 1
   }

   # Obtain subsampling directions

   if (ndk[3] > 0) {
      V = SdwDir(x0,ndk[3])           #####3ojoo!
      Vv = cbind(Vv,V)
   }

   # Undo standardization transformation

   Vv = solve(Rr,Vv)
   uaux = colSums(Vv^2)
   Vv = Vv %*% diag(1/sqrt(uaux))

   return(Vv)
}

################### Function MaxKur

MaxKur <- function(x,maxit = 30) {

   #
   # Compute direction maximizing the kurtosis of the projections
   # Observations passed to the routine should be standardized
   # Uses Newton's method and an augmented Lagrangian merit function
   # To be used as a subroutine of KurNwm
   #
   # Inputs:  x, observations (by rows)
   #          maxit, a limit to the number of Newton iterations
   # Outputs: v, max. kurtosis direction (local optimizer)
   #          f, max. kurtosis value
   #

   # Initialization

   ## Tolerances

   maxitini = 10
   mxitbl = 20
   tol = 1.0e-4
   tol1 = 1.0e-7
   tol2 = 1.0e-2
   tol3 = 1.0e-12
   beta = 1.0e-4
   rho0 = 0.1

   dx = dim(x)
   n = dx[1]
   p = dx[2]

   ep = matrix(1,p,1)

   ## Initial estimate of the direction

   uv = matrix(apply(x*x,1,sum),n,1)
   uw = 1/(tol3 + sqrt(uv))
   uu = x*(uw %*% t(ep))
   Su = cov(uu)
   Sue = eigen(Su)
   V = Sue$vectors
   D = Sue$values

   for (i in 1:p) {
      if (i == 1) r = ValKur(x,V[,i]) else r = cbind(r,ValKur(x,V[,i]))
   }
   ik = which.max(r)
   v = r[ik]
   a = V[,ik[1]]

   itini = 1
   difa = 1
   while ((itini <= maxitini)&&(difa > tol2)) {
      z = x %*% a
      zaux = z^2
      xaux = t(x)*(ep %*% t(zaux))
      H = 12*xaux %*% x
      He = eigen(H)
      V = He$vectors
      E = He$values
      iv = which.max(E)
      vv = E[iv]
      aa = V[,iv[1]]
      difa = norm(as.matrix(a - aa),type = "F")
      a = aa
      itini = itini + 1
   }
        vini=a
   ## Values at iteration 0 for the optimization algorithm

   z = x %*% a
   sk = sum(z^4)
   lam = 2*sk
   f = sk
   g = t(4*t(z^3) %*% x)
   zaux = z^2
   xaux = t(x)*(ep %*% t(zaux))
   H = 12*xaux %*% x

   al = 0
   it = 0
   diff = 1
   rho = rho0
   clkr = 0

   cc = 0

   # Newton method starting from the initial direction

   while (1) {

      ## Check termination conditions

      gl = g - 2*c(lam)*a

      A = 2*t(a)
      aqr = qr(a)
      Q = qr.Q(aqr, complete = T)
      Z = Q[,(2:p)]

      crit = norm(gl,type = "F") + abs(cc)
      if ((crit <= tol)||(it >= maxit)) break

      ## Compute search direction

      Hl = H - 2*lam[1]*diag(p)
      Hr = t(Z) %*% Hl %*% Z
      Hre = eigen(Hr)
      V = Hre$vectors
      E = Hre$values
      if (length(E) > 1) {
         Es1 = pmin(-abs(E),-tol)
         Hs = V %*% diag(Es1) %*% t(V)
      } else {
         Es1 = min(-abs(E),-tol)
         Hs = Es1 * V %*% t(V)
      }
      py = - cc/(A %*% t(A))
      rhs = t(Z) %*% (g + H %*% t(A) %*% py)
      pz = -solve(Hs,rhs)

      pp = Z %*% pz + t(A) %*% py

      dlam = (t(a) %*% (gl + H %*% pp))/(2*t(a) %*% a)

      ## Adjust penalty parameter

      f0d = t(gl) %*% pp - 2*rho*cc*t(a) %*% pp - dlam*cc
      crit1 = beta*norm(pp,type = "F")^2

      if (f0d < crit1) {
         rho1 = 2*(crit1 - f0d)/(tol3 + cc^2)
         rho = max(c(2*rho1,1.5*rho,rho0))
         f = sk - lam*cc - 0.5*rho*cc^2
         f0d = t(gl) %*% pp - 2*rho*cc*t(a) %*% pp - dlam*cc
         clkr = 0
      } else if ((f0d > 1000*crit1)&&(rho > rho0)) {
         rho1 = 2*(crit1 - t(gl) %*% pp + dlam*cc)/(tol3 + cc^2)
         if ((clkr == 4)&&(rho > 2*rho1)) {
            rho = 0.5*rho
            f = sk - lam*cc - 0.5*rho*cc^2
            f0d = t(gl) %*% pp - 2*rho*cc*t(a) %*% pp - dlam*cc
            clkr = 0
         } else clkr = clkr + 1
      }

      if (abs(f0d/(norm(g-2*rho*a*c(cc),type = "F")+abs(cc))) < tol1) break

      ## Line search

      al = 1
      itbl = 0
      while (itbl < mxitbl) {
         aa = a + al*pp
         lama = lam + al*dlam
         zz = x %*% aa
         cc = t(aa) %*% aa - 1
         sk = sum(zz^4)
         ff = sk - lama*cc - 0.5*rho*cc^2
         if (ff > f + 0.0001*al*f0d) break
         al = al/2
         itbl = itbl + 1
      }

      if (itbl >= mxitbl) break

      ## Update values for the next iteration

      a = aa
      lam = lama
      z = zz

      nmd2 = t(a) %*% a
      cc = nmd2 - 1
      f = sk - lam*cc - 0.5*rho*cc^2
      g = t(4*t(z^3) %*% x)
      zaux = z^2
      xaux = t(x)*(ep %*% t(zaux))
      H = 12*xaux %*% x

      it = it + 1
   }

   # Values to be returned

   a = matrix(c(a),p,1)
   v = a/norm(a,type = "F")
   xa = x %*% v
   f = sum(xa^4)/n
   xain = x %*% vini
   fini = sum(xain^4)/n
   rval = list(v = v, f = f, fini=fini,vini=vini)
   return(rval)
}

################### Function MinKur

MinKur <- function(x,maxit = 30) {

   #
   # Compute direction minimizing the kurtosis of the projections
   # Observations passed to the routine should be standardized
   # Uses Newton's method and an augmented Lagrangian merit function
   # To be used as a subroutine of KurNwm
   #
   # Inputs:  x, observations (by rows)
   #          maxit, a limit to the number of Newton iterations
   # Outputs: v, max. kurtosis direction (local optimizer)
   #          f, max. kurtosis value
   #

   # Initialization                        ########## Recheck all

   ## Tolerances

   maxitini = 10
   mxitbl = 20
   tol = 1.0e-4
   tol1 = 1.0e-7
   tol2 = 1.0e-2
   tol3 = 1.0e-12
   beta = 1.0e-4
   rho0 = 0.1

   dx = dim(x)
   n = dx[1]
   p = dx[2]

   ep = matrix(1,p,1)

   ## Initial estimate of the direction

   uv = matrix(apply(x*x,1,sum),n,1)
   uw = 1/(tol3 + sqrt(uv))
   uu = x*(uw %*% t(ep))

   Su = cov(uu)
   Sue = eigen(Su)
   V = Sue$vectors
   D = Sue$values

   for (i in 1:p) {
      if (i == 1) r = ValKur(x,V[,i]) else r = cbind(r,ValKur(x,V[,i]))
   }
   ik = which.min(r)
   v = r[ik]
   a = V[,ik[1]]

   itini = 1
   difa = 1
   while ((itini <= maxitini)&&(difa > tol2)) {
      z = x %*% a
      zaux = z^2
      xaux = t(x)*(ep %*% t(zaux))
      H = 12*xaux %*% x
      He = eigen(H)
      V = He$vectors
      E = He$values
      iv = which.min(E)
      vv = E[iv]
      aa = V[,iv[1]]
      difa = norm(as.matrix(a - aa),type = "F")
      a = aa
      itini = itini + 1
   }
   aini=a

   ## Values at iteration 0 for the optimization algorithm

   z = x %*% a
   sk = sum(z^4)
   lam = 2*sk
   f = sk
   g = t(4*t(z^3) %*% x)
   zaux = z^2
   xaux = t(x)*(ep %*% t(zaux))
   H = 12*xaux %*% x

   al = 0
   it = 0
   diff = 1
   rho = rho0
   clkr = 0

   cc = 0

   # Newton method starting from the initial direction

   while (1) {

   ## Check termination conditions

      gl = g - 2*c(lam)*a

      A = 2*t(a)
      aqr = qr(a)
      Q = qr.Q(aqr, complete = T)
      Z = Q[,(2:p)]

      crit = norm(gl,type = "F") + abs(cc)
      if ((crit <= tol)||(it >= maxit)) break

   ## Compute search direction

      Hl = H - 2*lam[1]*diag(p)
      Hr = t(Z) %*% Hl %*% Z
      Hre = eigen(Hr)
      V = Hre$vectors
      E = Hre$values
      if (length(E) > 1) {
         Es2 = pmax(abs(E),tol)
         Hs = V %*% diag(Es2) %*% t(V)
      } else {
         Es2 = max(abs(E),tol)
         Hs = Es2 * V %*% t(V)
      }
      py = - cc/(A %*% t(A))
      rhs = t(Z) %*% (g + H %*% t(A) %*% py)
      pz = -solve(Hs,rhs)
      pp = Z %*% pz + t(A) %*% py

      dlam = (t(a) %*% (gl + H %*% pp))/(2*t(a) %*% a)

   ## Adjust penalty parameter

      f0d = t(gl) %*% pp + 2*rho*cc*t(a) %*% pp - dlam*cc
      crit1 = beta*norm(pp,type = "F")^2

      if (f0d < crit1) {
         rho1 = 2*(crit1 + f0d)/(tol3 + cc^2)
         rho = max(c(2*rho1,1.5*rho,rho0))
         f = sk - lam*cc + 0.5*rho*cc^2
         f0d = t(gl) %*% pp + 2*rho*cc*t(a) %*% pp - dlam*cc
         clkr = 0
      } else if ((f0d < -1000*crit1)&&(rho > rho0)) {
         rho1 = 2*(crit1 - t(gl) %*% pp + dlam*cc)/(tol3 + cc^2)
         if ((clkr == 4)&&(rho > 2*rho1)) {
            rho = 0.5*rho
            f = sk - lam*cc + 0.5*rho*cc^2
            f0d = t(gl) %*% pp + 2*rho*cc*t(a) %*% pp - dlam*cc
            clkr = 0
         } else clkr = clkr + 1
      }

      if (abs(f0d/(norm(g-2*rho*a*c(cc),type = "F")+abs(cc))) < tol1) break

   ## Line search

      al = 1
      itbl = 0
      while (itbl < mxitbl) {
         aa = a + al*pp
         lama = lam + al*dlam
         zz = x %*% aa
         cc = t(aa) %*% aa - 1
         sk = sum(zz^4)
         ff = sk - lama*cc + 0.5*rho*cc^2
         if (ff < f + 0.0001*al*f0d) break
         al = al/2
         itbl = itbl + 1
      }

      if (itbl >= mxitbl) break

   ## Update values for the next iteration

      a = aa
      lam = lama
      z = zz

      nmd2 = t(a) %*% a
      cc = nmd2 - 1
      f = sk - lam*cc + 0.5*rho*cc^2
      g = t(4*t(z^3) %*% x)
      zaux = z^2
      xaux = t(x)*(ep %*% t(zaux))
      H = 12*xaux %*% x

      it = it + 1
   }

# Values to be returned

   a = matrix(c(a),p,1)
   v = a/norm(a,type = "F")
   xa = x %*% v
   f = sum(xa^4)/n
   vini=aini; xain=x%*%vini
   fini=mean(xain^4)
   rval = list(v = v, f = f, vini=vini,fini=fini)
   return(rval)
}

################### Function SdwDir

SdwDir <- function(x,nd) {

   #
   # Directions computed in a manner similar to Stahel-Donoho subsampling
   # Pairs of observations are chosen, and weights are assigned based on
   # the projections, by grouping the projected observations
   #
   #  Inputs:    x, observations by rows
   #             nd, number of directions to generate
   #
   #  Outputs:   V, directions generated by the procedure
   #

   use_rnd = 1

   # if (use_rnd == 0) {
   #    k_rnd = 1
   #    load("lstunif.RData")    # Load data defined as a variable (list) lstunif
   #    n_rnd = 100
   # }

   sx = dim(x)
   n = sx[1]
   p = sx[2]

   ep = matrix(1,p,1)

   ep5 = 1.0e-07
   mdi = floor((n+1)/2)
   k1 = p-1+floor((n+1)/2)
   k2 = p-1+floor((n+2)/2)
   fct = qnorm((k1/n + 1)/2)

   ## Number of groups to consider in the projections for any weight set

   s_grp = 2*p
   n_grp = floor(n/s_grp)

   # Generating the directions

   ndir = 0
   V = matrix(0,p,0)

   while (ndir < nd) {

      ## Sampling to obtain initial projections and weights

      # if (use_rnd == 1) {
         nn1 = min(1 + floor(runif(1)*n),n)
      # } else {
      #    nn1 = min(1 + floor(lstunif[k_rnd]*n),n)
      #    k_rnd = k_rnd + 1
      #    if (k_rnd > n_rnd) k_rnd = 1
      # }

      nn2 = nn1

      while ((nn2 == nn1)||(norm(as.matrix(x[nn1,1:p] - x[nn2,1:p]),type = "F") < ep5)) {
         # if (use_rnd == 1) {
            nn2 = min(1 + floor(runif(1)*n),n)
         # } else {
         #    nn2 = min(1 + floor(lstunif[k_rnd]*n),n)
         #    k_rnd = k_rnd + 1
         #    if (k_rnd > n_rnd) k_rnd = 1
         # }
      }
      dd = as.matrix(x[nn1,1:p] - x[nn2,1:p])
      dd = dd/norm(dd,type = "F")
      pr = x %*% dd
      prs = sort.int(pr, index.return = T)
      ww = prs$x
      lbl = prs$ix

      ## Stahel-Donoho directions

      j = 1
      while ((j <= n_grp)&&(ndir < nd)) {
         # if (use_rnd == 1) {
            auxs = sort.int(runif(s_grp), index.return = T)
         # } else {
         #    auxs = sort.int(lstunif[k_rnd:(k_rnd+s_grp-1)], index.return = T)
         #    k_rnd = k_rnd + s_grp
         #    if (k_rnd > n_rnd) k_rnd = 1
         # }
         bb = auxs$ix
         nn = bb[1:p]
         idr = lbl[(1+(j-1)*s_grp):(j*s_grp)]
         y1 = x[idr[nn],]
         qry1 = qr(y1)
         if (qry1$rank >= p) {
            dd = solve(qry1,ep)
            if (norm(dd,type = "F") > ep5) dd = dd/norm(dd,type = "F")
         }
         if (dim(V)[2] > 0) V = cbind(V,dd) else V = dd
         ndir = ndir + 1
         j = j + 1
      }
   }

   ## Return values

   return(V)
}

################### Function ValKur

ValKur <- function(x,d,km = 4) {

   #
   # Cluster identification from projections onto directions
   # maximizing and minimizing the kurtosis coefficient of
   # the data
   #
   # Evaluate the moment coefficient of order k
   # for the univariate projection of multivariate data
   #
   # Inputs:  x, observations (by rows)
   #          d, projection direction
   #          km, order of the moment (km = 4 by default)
   # Output:  mcv, value of the moment coefficient for the univariate data

   dimx = dim(x)
   n = dimx[1]

   t = x %*% d
   tm = mean(t)
   tt = abs(t - tm)
   vr = sum(tt^2)/(n-1)
   kr = sum(tt^km)/n
   mcv = kr/(vr^(km/2))

   ## Return values

   return(mcv)
}

