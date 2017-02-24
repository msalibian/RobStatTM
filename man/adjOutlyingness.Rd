\name{adjOutlyingness}
\alias{adjOutlyingness}
\title{Compute (Skewness-adjusted) Multivariate Outlyingness}
\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\description{
  For an \eqn{n \times p}{n * p} data matrix (or data frame) \code{x},
  compute the \dQuote{\emph{outlyingness}} of all \eqn{n} observations.
  Outlyingness here is a generalization of the Donoho-Stahel
  outlyingness measure, where skewness is taken into account via the
  medcouple, \code{\link{mc}()}.
}
\usage{
adjOutlyingness(x, ndir = 250, clower = 4, cupper = 3,
                alpha.cutoff = 0.75, coef = 1.5,
                qr.tol = 1e-12, keep.tol = 1e-12,
                only.outlyingness = FALSE, maxit.mult = max(100, p))
}
\arguments{
  \item{x}{a numeric \code{\link{matrix}} or \code{\link{data.frame}}.}
  \item{ndir}{positive integer specifying the number of directions that
    should be searched.}
  \item{clower, cupper}{the constant to be used for the lower and upper
    tails, in order to transform the data towards symmetry.  You can set
    \code{clower = 0, cupper = 0} to get the \emph{non}-adjusted,
    i.e., classical (\dQuote{central} or \dQuote{symmetric})
    outlyingness.  In that case, \code{\link{mc}()} is not used.}
  \item{alpha.cutoff}{number in (0,1) specifying the quantiles
    \eqn{(\alpha, 1-\alpha)} which determine the \dQuote{outlier}
    cutoff.  The default, using quartiles, corresponds to the definition
    of the medcouple (\code{\link{mc}}), but there is no stringent
    reason for using the same alpha for the outlier cutoff.}
  \item{coef}{positive number specifying the factor with which the
    interquartile range (\code{\link{IQR}}) is multiplied to determine
    \sQuote{boxplot hinges}-like upper and lower bounds.}
  \item{qr.tol}{positive tolerance to be used for \code{\link{qr}} and
    \code{\link{solve.qr}} for determining the \code{ndir} directions,
    each determined by a random sample of \eqn{p} (out of \eqn{n})
    observations.}
  \item{keep.tol}{positive tolerance to determine which of the sample
    direction should be kept, namely only those for which
    \eqn{\|x\| \cdot \|B\|}{||x|| * ||B||} is larger than \code{keep.tol}.}
  \item{only.outlyingness}{logical indicating if the final outlier
    determination should be skipped.  In that case, a vector is
    returned, see \sQuote{Value:} below.}
  \item{maxit.mult}{integer factor; \code{maxit <- maxit.mult * ndir}
    will determine the maximal number of direction searching
    iterations.  May need to be increased for higher dimensional data.}
}
\note{
 The result is \emph{random} as it depends on the sample of
 \code{ndir} directions chosen.

 Till Aug/Oct. 2014, the default values for \code{clower} and \code{cupper} were
 accidentally reversed, and the signs inside \code{exp(.)} where swapped
 in the (now corrected) two expressions \preformatted{
 tup <- Q3 + coef * IQR * exp(.... + clower * tmc * (tmc < 0))
 tlo <- Q1 - coef * IQR * exp(.... - cupper * tmc * (tmc < 0))
}
 already in the code from Antwerpen (\file{mcrsoft/adjoutlingness.R}),
 contrary to the published reference.

 Further, the original algorithm had not been scale-equivariant in the
 direction construction, which has been amended in 2014-10 as well.

 The results, including diagnosed outliers, therefore have changed,
 typically slightly, since \pkg{robustbase} version 0.92-0.
}
\details{
  \bold{FIXME}:  Details in the comment of the Matlab code;
  also in the reference(s).
  %% SEE /u/maechler/R/MM/STATISTICS/robust/MC/mcmatl/adjoutlyingness.m
  %% ---- which has notes about input/output etc of the corresponding
  %%      Matlab code

  The method as described can be useful as preprocessing in
  FASTICA (\url{http://www.cis.hut.fi/projects/ica/fastica/};
  see also the \R package \CRANpkg{fastICA}.
}
\value{
  If \code{only.outlyingness} is true, a vector \code{adjout},
  otherwise, as by default, a list with components
  \item{adjout}{numeric of \code{length(n)} giving the adjusted
    outlyingness of each observation.}
  \item{cutoff}{cutoff for \dQuote{outlier} with respect to the adjusted
    outlyingnesses, and depending on \code{alpha.cutoff}.}
  \item{nonOut}{logical of \code{length(n)}, \code{TRUE} when the
    corresponding observation is \bold{non}-outlying with respect to the
    cutoff and the adjusted outlyingnesses.}
}
\references{
  Brys, G., Hubert, M., and Rousseeuw, P.J. (2005)
  A Robustification of Independent Component Analysis;
  \emph{Journal of Chemometrics}, \bold{19}, 1--12.

  For the up-to-date reference, please consult
  \url{http://wis.kuleuven.be/stat/robust}
}
\author{Guy Brys; help page and improvements by Martin Maechler}
\seealso{the adjusted boxplot, \code{\link{adjbox}} and the medcouple,
  \code{\link{mc}}.
}
\examples{
## An Example with bad condition number and "border case" outliers

dim(longley)
set.seed(1) ## result is random!
ao1 <- adjOutlyingness(longley)
## which are outlying ?
which(!ao1$nonOut) ## one: "1948" - for this seed! (often: none)
stopifnot(all(ao1$nonOut[-2]))

## An Example with outliers :

dim(hbk)
set.seed(1)
ao.hbk <- adjOutlyingness(hbk)
str(ao.hbk)
hist(ao.hbk $adjout)## really two groups
table(ao.hbk$nonOut)## 14 outliers, 61 non-outliers:
## outliers are :
which(! ao.hbk$nonOut) # 1 .. 14   --- but not for all random seeds!

## here, they are the same as found by (much faster) MCD:
cc <- covMcd(hbk)
stopifnot(all(cc$mcd.wt == ao.hbk$nonOut))

## This is revealing (about 1--2 cases, where outliers are *not* == 1:14
##  but needs almost 1 [sec] per call:
if(interactive()) {
  for(i in 1:30) {
    print(system.time(ao.hbk <- adjOutlyingness(hbk)))
    if(!identical(iout <- which(!ao.hbk$nonOut), 1:14)) {
	 cat("Outliers:\n"); print(iout)
    }
  }
}

## "Central" outlyingness: *not* calling mc()  anymore, since 2014-12-11:
trace(mc)
out <- capture.output(
  oo <- adjOutlyingness(hbk, clower=0, cupper=0)
)
untrace(mc)
stopifnot(length(out) == 0)
}
\keyword{robust}
\keyword{multivariate}
