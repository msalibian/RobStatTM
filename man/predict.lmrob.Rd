\name{predict.lmrob}
\alias{predict.lmrob}
\title{Predict method for Robust Linear Model ("lmrob") Fits}
\description{
  Predicted values based on robust linear model object.
}
\usage{
\method{predict}{lmrob}(object, newdata, se.fit = FALSE,
       scale = NULL, df = NULL,
       interval = c("none", "confidence", "prediction"), level = 0.95,
       type = c("response", "terms"), terms = NULL,
       na.action = na.pass, pred.var = res.var/weights, weights = 1, ...)
}
\arguments{
%% the following is +- copy-pasted from  predict.lm.Rd:
  \item{object}{object of class inheriting from \code{"lmrob"}}
  \item{newdata}{an optional data frame in which to look for variables with
    which to predict.  If omitted, the fitted values are used.}
  \item{se.fit}{a switch indicating if standard errors are required.}
  \item{scale}{scale parameter for std.err. calculation}
  \item{df}{degrees of freedom for scale}
  \item{interval}{type of interval calculation.}
  \item{level}{tolerance/confidence level}
  \item{type}{Type of prediction (response or model term).}
  \item{terms}{if \code{type="terms"}, which terms (default is all terms)}
  \item{na.action}{function determining what should be done with missing
    values in \code{newdata}.  The default is to predict \code{NA}.}
  \item{pred.var}{the variance(s) for future observations to be assumed
    for prediction intervals.  See \sQuote{Details}.}
  \item{weights}{variance weights for prediction. This can be a numeric
    vector or a one-sided model formula. In the latter case, it is
    interpreted as an expression evaluated in \code{newdata}}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  Note that this \code{lmrob} method for \code{\link{predict}} is
  closely modeled after the method for \code{lm()},
  \code{\link{predict.lm}}, maybe see there for caveats with missing
  value treatment.
  %% Also lifted from  predict.lm.Rd :
  The prediction intervals are for a single observation at each case in
  \code{newdata} (or by default, the data used for the fit) with error
  variance(s) \code{pred.var}.  This can be a multiple of \code{res.var},
  the estimated value of \eqn{\sigma^2}: the default is to assume that
  future observations have the same error variance as those
  used for fitting.  If \code{weights} is supplied, the inverse of this
  is used as a scale factor.  For a weighted fit, if the prediction
  is for the original data frame, \code{weights} defaults to the weights
  used for the  model fit, with a warning since it might not be the
  intended result.  If the fit was weighted and \code{newdata} is given, the
  default is to assume constant prediction variance, with a warning.
}
\value{
%% the following is +- copy-pasted from  predict.lm.Rd:
  \code{predict.lmrob} produces a vector of predictions or a matrix of
  predictions and bounds with column names \code{fit}, \code{lwr}, and
  \code{upr} if \code{interval} is set.  If \code{se.fit} is
  \code{TRUE}, a list with the following components is returned:
  \item{fit}{vector or matrix as above}
  \item{se.fit}{standard error of predicted means}
  \item{residual.scale}{residual standard deviations}
  \item{df}{degrees of freedom for residual}
}
% \references{

% }
\author{Andreas Ruckstuhl}
\seealso{
  \code{\link{lmrob}} and the (non-robust) traditional
  \code{\link{predict.lm}} method.
}
\examples{
## Predictions --- artificial example -- closely following  example(predict.lm)

set.seed(5)
n <- length(x <- sort(c(round(rnorm(25), 1), 20)))
y <- x + rnorm(n)
iO <- c(sample(n-1, 3), n)
y[iO] <- y[iO] + 10*rcauchy(iO)

p.ex <- function(...) {
  plot(y ~ x, ...); abline(0,1, col="sky blue")
  points(y ~ x, subset=iO, col="red", pch=2)
  abline(lm   (y ~ x), col = "gray40")
  abline(lmrob(y ~ x), col = "forest green")
  legend("topleft", c("true", "Least Squares", "robust"),
         col = c("sky blue", "gray40", "forest green"), lwd=1.5, bty="n")
}
p.ex()

fm <- lmrob(y ~ x)
predict(fm)
new <- data.frame(x = seq(-3, 10, 0.25))
str(predict(fm, new, se.fit = TRUE))
pred.w.plim <- predict(fm, new, interval = "prediction")
pred.w.clim <- predict(fm, new, interval = "confidence")
pmat <- cbind(pred.w.clim, pred.w.plim[,-1])

matlines(new$x, pmat, lty = c(1,2,2,3,3))# add to first plot
## show zoom-in region :
rect(xleft = -3, ybottom = -20, xright = 10, ytop = 40,
     lty = 3, border="orange4")

## now zoom in :
p.ex(xlim = c(-3,10), ylim = c(-20, 40))
matlines(new$x, pmat, lty = c(1,2,2,3,3))
box(lty = 3, col="orange4", lwd=3)
legend("bottom", c("fit", "lwr CI", "upr CI", "lwr Pred.I", "upr Pred.I"),
       col = 1:5, lty=c(1,2,2,3,3), bty="n")

## Prediction intervals, special cases
##  The first three of these throw warnings
w <- 1 + x^2
fit <- lmrob(y ~ x)
wfit <- lmrob(y ~ x, weights = w)
predict(fit,       interval = "prediction")
predict(wfit,      interval = "prediction")
predict(wfit, new, interval = "prediction")
predict(wfit, new, interval = "prediction", weights = (new$x)^2) -> p.w2
p.w2
stopifnot(identical(p.w2, ## the same as using formula:
     predict(wfit, new, interval = "prediction", weights = ~x^2)))
}
\keyword{robust}
\keyword{regression}
