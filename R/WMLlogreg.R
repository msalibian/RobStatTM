#' Weighted likelihood estimator for the logistic model
#'
#' This function computes a weighted likelihood estimator for the logistic model, where
#' the weights penalize high leverage observations. In this version the weights are zero or one.
#'
#' @export WMLlogreg logregWML
#' @aliases WMLlogreg logregWML
#' @rdname WMLlogreg
#'
#' @param x0 p x n matrix of explanatory variables, p is the number of explanatory variables, n is the number of observations
#' @param y response vector
#' @param intercept 1 or 0 indicating if an intercept is included or or not
#'
#' @return A list with the following components:
#' \item{coefficients}{vector of regression coefficients}
#' \item{standard.deviation}{standard deviations of the regression coefficient estimators}
#' \item{fitted.values}{vector with the probabilities of success}
#' \item{residual.deviances}{residual deviances}
#' \item{cov}{covariance matrix of the regression estimates}
#' \item{objective}{value of the objective function at the minimum}
#' \item{xweights}{vector of zeros and ones used to compute the weighted maimum likelihood estimator}
#'
#' @author Victor Yohai
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(skin)
#' Xskin <- as.matrix( skin[, 1:2] )
#' yskin <- skin$vasoconst
#' skinWML <- logregWML(Xskin, yskin, intercept=1)
#' skinWML$coeff
#' skinWML$standard.deviation
#'
logregWML <- WMLlogreg <- function (x0, y, intercept = 1)
  {
  ttx=colnames(x0)
  if (!is.numeric(y))
        y <- as.numeric(y)
  if (!is.null(dim(y))) {
  if (ncol(y) != 1)
            stop("y is not onedimensional")
        y <- as.vector(y)
    }
    n <- length(y)
  #  if (is.data.frame(x0)) {
  #      x0 <- data.matrix(x0)
  #  }
  #else if (!is.matrix(x0)) {
  #    x0 <- matrix(x0, length(x0), 1, dimnames = list(names(x0),
  #       deparse(substitute(x0))))
  # }
   x0<-as.matrix(x0)
   if (nrow(x0) != n)
        stop("Number of observations in x and y not equal")
   na.x <- !is.finite(rowSums(x0))
   na.y <- !is.finite(y)
   ok <- !(na.x | na.y)
   if (!all(ok)) {
   x0 <- x0[ok, , drop = FALSE]
   y <- y[ok]
    }
   n <- nrow(x0)
   if (n == 0)
   stop("All observations have missing values!")
   p<-ncol(x0)
   family <- binomial()
   p<-dim(x0)[2]
   zz<-rep(0,p)
   for (i in 1:p)
   {zz[i]<-sum((x0[,i]==0)|(x0[,i]==1))}
   tt<-which(zz!=n)
   p1<-length(tt)
   x0=as.matrix(x0,n,p)
   x00<-x0[,tt]
   if(p1==1)
   {rdx<-abs(x00-median(x00))/mad(x00)
   wrd<-(rdx<=qnorm(.9875))}
   if(p1>1)
   {
   mcdx<-robustbase::covMcd(x00,alpha=.75)
   rdx<-mahalanobis(x00,center=mcdx$center,cov=mcdx$cov)
   vc<-qchisq(0.975,p)
   wrd<-(rdx<=vc)}
   if(p1==0)
   {wrd=1:n}

   x000=x0[wrd,]
   colnames(x000)=ttx
   y00=y[wrd]

   if (intercept==1)
   {
   out=glm(y00~x000, family = family)
    x <- cbind(Intercept = 1, x0)}
   if (intercept==0)
   {out<-glm(y00~x000-1, family = family)
   x<-x0}
   out1<-summary(out)
   fitted.linear<-(x%*%out$coeff)
   fitted.linear<-pmin(fitted.linear,1e2)
   fitted.values<-exp(fitted.linear)/(1+exp(fitted.linear))
   residual.deviances<-sign(y-fitted.values)*sqrt(-2*(y*log(fitted.values)+(1-y)*log(1-fitted.values)))
   result<-list(xweights=wrd, coefficients=out$coeff,standard.deviation=sqrt(diag(out1$cov.uns)),fitted.values=t(fitted.values),cov=out1$cov.unscaled, residual.deviances= t(residual.deviances))
   result
 }




