\name{summary.glmrob}
\alias{summary.glmrob}
\alias{vcov.glmrob}
\alias{print.summary.glmrob}
\title{Summarizing Robust Fits of Generalized Linear Models}
\description{
  The \code{summary} method for class \code{"\link{glmrob}"}
  summarizes robust fits of (currently only discrete) generalized linear
  models.
}
\usage{
\method{summary}{glmrob}(object, correlation = FALSE, symbolic.cor = FALSE, \dots)
\method{vcov}{glmrob}(object, \dots)

\method{print}{summary.glmrob}(x, digits = max(3, getOption("digits") - 3),
      symbolic.cor = x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"glmrob"}, usually, a result of
    a call to \code{\link{glmrob}}.}
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned and printed.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in
    a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
  \item{\dots}{further arguments passed to or from other methods.}
  \item{x}{an object of class \code{"summary.glrob"}.}
  \item{digits}{the number of digits to use for printing.}
  \item{signif.stars}{logical indicating if the P-values should be
    visualized by so called \dQuote{significance stars}.}
}

\details{
  \code{\link{summary.glmrob}} returns an object of \code{\link{class}}
  \code{"summary.glmrob"}.

  Its \code{\link{print}()} method tries to be smart about formatting the
  coefficients, standard errors, etc, and gives
  \dQuote{significance stars} if \code{signif.stars} is \code{TRUE}
  (as per default when \code{\link{options}} where not changed).
}

\value{
  The function \code{\link{summary.glmrob}} computes and returns a list
  of summary statistics of the robustly fitted linear model given in
  \code{object}. The following elements are in the list:
  \item{...}{ FIXME }
  %% MM thinks that summary should possibly just *add* to the original
  %% object.
}

%%\references{ ~put references to the literature/web site here ~ }
\author{Andreas Ruckstuhl}

%%\note{ ~~further notes~~ }
%%
\seealso{\code{\link{glmrob}}; the generic \code{\link{summary}} and
  also \code{\link[stats]{summary.glm}}.
}
\examples{
data(epilepsy)
Rmod <- glmrob(Ysum ~ Age10 + Base4*Trt, family = poisson,
               data = epilepsy, method= "Mqle")
ss <- summary(Rmod)
ss ## calls print.summary.glmrob()
str(ss) ## internal STRucture of summary object
}
\keyword{regression}
\keyword{nonlinear}
\keyword{robust}
