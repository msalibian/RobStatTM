\name{plot.lts}
\alias{plot.lts}
\alias{ltsPlot}
\title{Robust LTS Regression Diagnostic Plots}
\description{
  Four plots (selectable by \code{which}) are currently provided:
  \enumerate{
    \item a plot of the standardized residuals versus their index,
    \item a plot of the standardized residuals versus fitted values,
    \item a Normal Q-Q plot of the standardized residuals, and
    \item a regression diagnostic plot (standardized residuals versus
      robust distances of the predictor variables).
  }
}
\usage{
\method{plot}{lts}(x, which = c("all","rqq","rindex","rfit","rdiag"),
     classic=FALSE, ask=(which=="all" && dev.interactive()), id.n, \dots)
%    ltsPlot(x, which = c("all","rqq","rindex","rfit","rdiag"), classic=FALSE, ask=(which=="all" && dev.interactive()), id.n, \dots)
}
\arguments{
  \item{x}{a \code{lts} object, typically result of \code{ltsReg}.}
  \item{which}{string indicating which plot to show.  See the
    \emph{Details} section for a description of the options.  Defaults
    to \code{"all"}.}.
  \item{classic}{whether to plot the classical distances too. Default is
    \code{FALSE}.}.
  \item{ask}{logical indicating if the user should be \emph{ask}ed
    before each plot, see \code{\link{par}(ask=.)}.  Defaults to
    \code{which == "all" && \link{dev.interactive}()}.
  }
  \item{id.n}{number of observations to be identified by a label starting
    with the most extreme.  Default is the number of identified outliers
    (can be different for the different plots - see Details).}
  \item{\dots}{other parameters to be passed through to plotting functions.}
}
\details{
  This function produces several plots based on the robust and classical
  regression estimates.
  Which of them to select is specified by the attribute  \code{which}.
  The possible options are:
  \describe{
    \item{\code{rqq}:}{Normal Q-Q plot of the standardized residuals;}
    \item{\code{rindex}:}{plot of the standardized residuals versus their
      index;}
    \item{\code{rfit}:}{plot of the standardized residuals versus fitted
      values;}
    \item{\code{rdiag}:}{regression diagnostic plot.}
  }
  The normal quantile plot produces a normal Q-Q plot of the
  standardized residuals.
  A line is drawn which passes through the first and third quantile.
  The \code{id.n} residuals with largest distances from this line are
  identified by labels (the observation number).  The default for
  \code{id.n} is the number of regression outliers (lts.wt==0).

  In the Index plot and in the Fitted values plot the standardized
  residuals are displayed against the observation number or the fitted
  value respectively.
  A horizontal dashed line is drawn at 0 and two solid horizontal lines are
  located at +2.5 and -2.5. The id.n residuals with largest absolute values are
  identified by labels (the observation number).  The default for id.n is the
  number regression outliers (lts.wt==0).

  The regression diagnostic plot, introduced by Rousseeuw and van
  Zomeren (1990), displays the standardized residuals versus robust
  distances. Following Rousseeuw and van Zomeren (1990), the
  horizontal dashed lines are located at +2.5 and -2.5  and the
  vertical line is located at the upper 0.975 percent point of the
  chi-squared distribution with p degrees of freedom. The id.n residuals
  with largest absolute values and/or largest robust Mahalanobis distances are
  identified by labels (the observation number). The default for id.n is
  the number of all outliers: regression outliers (lts.wt==0) + leverage
  (bad and good) points (RD > 0.975 percent point of the chi-squared
  distribution with p degrees of freedom).
}
%\value{}
\references{
  P. J. Rousseeuw and van Zomeren, B. C. (1990).
  Unmasking Multivariate Outliers and Leverage Points.
  \emph{Journal of the American Statistical Association} \bold{85}, 633--639.

  P. J. Rousseeuw and K. van Driessen (1999)
  A fast algorithm for the minimum covariance determinant estimator.
  \emph{Technometrics} \bold{41}, 212--223.
}
\seealso{
  \code{\link{covPlot}}
}
\examples{
data(hbk)
lts <- ltsReg(Y ~ ., data = hbk)
lts
plot(lts, which = "rqq") %% Improve: also show model 'call'
}
\keyword{hplot}
\keyword{multivariate}
