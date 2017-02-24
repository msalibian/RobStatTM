\name{splitFrame}
\alias{splitFrame}
\title{
  Split Continuous and Categorical Predictors
}
\description{
  Splits the design matrix into categorical and continuous
  predictors.  Categorical variables are variables that are factors or
  ordered factors.
}
\usage{
splitFrame(mf, x = model.matrix(mt, mf),
 	   type = c("f","fi", "fii"))
}
\arguments{
  \item{mf}{model frame (as returned by \code{\link{model.frame}}).}
  \item{x}{(optional) design matrix, defaulting to the derived
    \code{\link{model.matrix}}.}
  \item{type}{a character string specifying the split type (see details).}
}
\details{
  Which split type is used can be controlled with the setting
  \code{split.type} in \code{\link{lmrob.control}}.

  There are three split types. The only differences between the types
  are how interactions between categorical and continuous variables are
  handled. The extra types of splitting can be used to avoid
  \emph{Too many singular resamples} errors.

  Type \code{"f"}, the default, assigns only the intercept, categorical and
  interactions of categorical variables to \code{x1}.  Interactions of
  categorical and continuous variables are assigned to \code{x2}.

  Type \code{"fi"} assigns also interactions between categorical and
  continuous variables to \code{x1}.

  Type \code{"fii"} assigns not only interactions between categorical and
  continuous variables to \code{x1},  but also the (corresponding)
  continuous variables themselves.
}
\value{
  A list that includes the following components:
  \item{x1 }{design matrix containing only categorical variables}
  \item{x1.idx }{logical vectors of the variables considered
    categorical in the original design matrix}
  \item{x2 }{design matrix containing the continuous variables}
}
\references{
  Maronna, R. A., and Yohai, V. J. (2000).
  Robust regression with both continuous and categorical predictors.
  \emph{Journal of Statistical Planning and Inference} \bold{89}, 197--214.
}
\author{
  Manuel Koller
}
\seealso{
  \code{\link{lmrob.M.S}}
}
\examples{
data(education)
education <- within(education, Region <- factor(Region))

## no interactions -- same split for all types:
fm1 <- lm(Y ~ Region + X1 + X2 + X3, education)
splt <- splitFrame(fm1$model)
str(splt)

## with interactions:
fm2 <- lm(Y ~ Region:X1:X2 + X1*X2, education)
s1 <- splitFrame(fm2$model, type="f"  )
s2 <- splitFrame(fm2$model, type="fi" )
s3 <- splitFrame(fm2$model, type="fii")
cbind(s1$x1.idx,
      s2$x1.idx,
      s3$x1.idx)
rbind(p.x1 = c(ncol(s1$x1), ncol(s2$x1), ncol(s3$x1)),
      p.x2 = c(ncol(s1$x2), ncol(s2$x2), ncol(s3$x2)))
}
