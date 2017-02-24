\name{ltsReg}
\alias{ltsReg}
\alias{ltsReg.default}
\alias{ltsReg.formula}
\alias{print.lts}
\title{Least Trimmed Squares Robust (High Breakdown) Regression}
\concept{High breakdown point}
\description{
  Carries out least trimmed squares (LTS) robust (high breakdown point)
  regression.
}
\usage{
ltsReg(x, \dots)

\method{ltsReg}{formula}(formula, data, subset, weights, na.action,
       model = TRUE, x.ret = FALSE, y.ret = FALSE,
       contrasts = NULL, offset, \dots)

\method{ltsReg}{default}(x, y, intercept = TRUE, alpha = , nsamp = , adjust = ,
       mcd = TRUE, qr.out = FALSE, yname = NULL,
       seed = , trace = , use.correction = , wgtFUN = , control = rrcov.control(),
       \dots)
}
\arguments{
  \item{formula}{a \code{\link{formula}} of the form \code{y ~ x1 + x2 + ...}.}
  \item{data}{data frame from which variables specified in
    \code{formula} are to be taken.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{weights}{an optional vector of weights to be used
    in the fitting process. \bold{NOT USED YET}.
    %%% If specified, weighted least squares is used
    %%% with weights \code{weights} (that is, minimizing \code{sum(w*e^2)});
    %%% otherwise ordinary least squares is used.
    }
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \dQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{model, x.ret, y.ret}{\code{\link{logical}}s indicating if the
    model frame, the model matrix and the response are to be returned,
    respectively.}

  \item{contrasts}{an optional list.  See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.}
  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.  An \code{\link{offset}} term can be included in the
    formula instead or as well, and if both are specified their sum is used.}

  \item{x}{a matrix or data frame containing the explanatory variables.}
  \item{y}{the response: a vector of length the number of rows of \code{x}.}.
  \item{intercept}{if true, a model with constant term will be
  estimated; otherwise no constant term will be included.  Default is
  \code{intercept = TRUE}  }
  \item{alpha}{the percentage (roughly) of squared residuals whose sum will be
    minimized, by default 0.5.  In general, \code{alpha} must between
    0.5 and 1.}
  \item{nsamp}{number of subsets used for initial estimates or
    \code{"best"} or \code{"exact"}.  Default is \code{nsamp = 500}.  For
    \code{nsamp="best"} exhaustive enumeration is done, as long as the
    number of trials does not exceed 5000.  For \code{"exact"},
    exhaustive enumeration will be attempted however many samples are needed.
    In this case a warning message will be displayed saying that the
    computation can take a very long time. }
  \item{adjust}{whether to perform intercept adjustment at each step.
    Since this can be time consuming, the default is \code{adjust = FALSE}.}
  \item{mcd}{whether to compute robust distances using Fast-MCD.}
  \item{qr.out}{whether to return the QR decomposition (see
    \code{\link{qr}}); defaults to false.}
  \item{yname}{the name of the dependent variable.  Default is \code{yname = NULL}}
  \item{seed}{initial seed for random generator, like
    \code{\link{.Random.seed}}, see \code{\link{rrcov.control}}.}
  \item{trace}{logical (or integer) indicating if intermediate results
    should be printed; defaults to \code{FALSE}; values \eqn{\ge 2}{>= 2}
    also produce print from the internal (Fortran) code.}
  \item{use.correction}{ whether to use finite sample correction factors.
    Default is \code{use.correction=TRUE}}
  \item{wgtFUN}{a character string or \code{\link{function}}, specifying
    how the weights for the reweighting step should be computed.
    Up to April 2013, the only option has been the original proposal in
    (1999), now specified by \code{wgtFUN = "01.original"} (or via \code{control}).}
  %% MM:  want control also for formula !?!?!?!?!!
  \item{control}{a list with estimation options - same as these provided
    in the function specification.  If the control object is supplied, the
    parameters from it will be used.  If parameters are passed also in the
    invocation statement, they will override the corresponding elements of
    the control object.}
  \item{\dots}{arguments passed to or from other methods.}
}
\details{
  The LTS regression method minimizes the sum of the \eqn{h} smallest
  squared residuals, where \eqn{h > n/2}, i.e. at least half the number of
  observations must be used.  The default value of \eqn{h} (when
  \code{alpha=1/2}) is roughly \eqn{n / 2}, more precisely,
  \code{(n+p+1) \%/\% 2} where \eqn{n} is the
  total number of observations, but by setting \code{alpha}, the user
  may choose higher values up to n, where
  \eqn{h = h(\alpha,n,p) =} \code{\link{h.alpha.n}(alpha,n,p)}.  The LTS
  estimate of the error scale is given by the minimum of the objective
  function multiplied by a consistency factor
  and a finite sample correction factor -- see Pison et al. (2002)
  for details.  The rescaling factors for the raw and final estimates are
  returned also in the vectors \code{raw.cnp2} and \code{cnp2} of
  length 2 respectively.  The finite sample corrections can be suppressed
  by setting \code{use.correction=FALSE}.  The computations are performed
  using the Fast LTS algorithm proposed by Rousseeuw and Van Driessen (1999).

  As always, the formula interface has an implied intercept term which can be
  removed either by \code{y ~ x - 1} or \code{y ~ 0 + x}.  See
  \code{\link{formula}} for more details.
}
\note{We strongly recommend using \code{\link{lmrob}()} instead of
  \code{ltsReg} (\emph{See also} below)!
}
\value{
  The function \code{ltsReg} returns an object of class \code{"lts"}.
  The \code{\link{summary}} method function is used to obtain (and
  print) a summary table of the results, and \code{\link[=ltsPlot]{plot}()}
  can be used to plot them, see the the specific help pages.

  The generic accessor functions \code{\link{coefficients}},
  \code{\link{fitted.values}} and \code{\link{residuals}}
  extract various useful features of the value returned by
  \code{ltsReg}.

  An object of class \code{lts} is a \code{\link{list}} containing at
  least the following components:
  \item{crit}{
    the value of the objective function of the LTS regression method,
    i.e., the sum of the \eqn{h} smallest squared raw residuals.
  }
  \item{coefficients}{
    vector of coefficient estimates (including the intercept by default when
    \code{intercept=TRUE}), obtained after reweighting.
  }
  \item{best}{
    the best subset found and used for computing the raw estimates, with
    \code{\link{length}(best) == quan = \link{h.alpha.n}(alpha,n,p)}.
  }
  \item{fitted.values}{vector like \code{y} containing the fitted values
    of the response after reweighting.}
  \item{residuals}{vector like \code{y} containing the residuals from
    the weighted least squares regression.}
  \item{scale}{scale estimate of the reweighted residuals.  }
  \item{alpha}{same as the input parameter \code{alpha}.}
  \item{quan}{the number \eqn{h} of observations which have determined
    the least trimmed squares estimator.}
  \item{intercept}{same as the input parameter \code{intercept}.}
  \item{cnp2}{a vector of length two containing the consistency
    correction factor and the finite sample correction factor of
    the final estimate of the error scale.}
  \item{raw.coefficients}{vector of raw coefficient estimates (including
  the intercept, when \code{intercept=TRUE}).}
  \item{raw.scale}{scale estimate of the raw residuals.}
  \item{raw.resid}{vector like \code{y} containing the raw residuals
    from the regression.}
  \item{raw.cnp2}{a vector of length two containing the consistency
    correction factor and the finite sample correction factor of the
    raw estimate of the error scale.}
  \item{lts.wt}{
    vector like y containing weights that can be used in a weighted
    least squares.  These weights are 1 for points with reasonably
    small residuals, and 0 for points with large residuals.
  }
  \item{raw.weights}{
    vector containing the raw weights based on the raw residuals and raw scale.
  }
  \item{method}{character string naming the method (Least Trimmed Squares).}
  \item{X}{the input data as a matrix (including intercept column if
    applicable).}
  \item{Y}{the response variable as a vector.}

}
\author{Valentin Todorov \email{valentin.todorov@chello.at}, based on
  work written for S-plus by Peter Rousseeuw and Katrien van Driessen
  from University of Antwerp.% no E-mails for spam-protection
}

\references{
  Peter J. Rousseeuw (1984), Least Median of Squares Regression.
  \emph{Journal of the American Statistical Association} \bold{79}, 871--881.

  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection.} Wiley.

  P. J. Rousseeuw and K. van Driessen (1999)
  A fast algorithm for the minimum covariance determinant estimator.
  \emph{Technometrics} \bold{41}, 212--223.

  Pison, G., Van Aelst, S., and Willems, G. (2002)
  Small Sample Corrections for LTS and MCD.
  \emph{Metrika} \bold{55}, 111-123.
}

\seealso{
  \code{\link{lmrob.S}()} provides a fast S estimator with similar
  breakdown point as \code{ltsReg()} but better efficiency.\cr
  For data analysis, rather use \code{\link{lmrob}} which is based on
  \code{\link{lmrob.S}}.

  \code{\link{covMcd}};
  \code{\link{summary.lts}} for summaries.

  The generic functions \code{\link{coef}}, \code{\link{residuals}},
  \code{\link{fitted}}.
}
\examples{
data(heart)
## Default method works with 'x'-matrix and y-var:
heart.x <- data.matrix(heart[, 1:2]) # the X-variables
heart.y <- heart[,"clength"]
ltsReg(heart.x, heart.y)

data(stackloss)
ltsReg(stack.loss ~ ., data = stackloss)
}
\keyword{robust}
\keyword{regression}
