\name{carrots}
\alias{carrots}
\docType{data}
\title{Insect Damages on Carrots}
\description{
  The damage carrots data set from Phelps (1982) was used by McCullagh
  and Nelder (1989) in order to illustrate diagnostic techniques because
  of the presence of an outlier. In a soil experiment trial with three
  blocks, eight levels of insecticide were applied and the carrots were
  tested for insect damage.
}
\usage{data(carrots)}

\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{success}{ integer giving the number of carrots with insect damage.}
    \item{total}{ integer giving the total number of carrots per
      experimental unit.}
    \item{logdose}{a numeric vector giving log(dose) values (eight
      different levels only).}
    \item{block}{factor with levels \code{B1} to \code{B3}}
  }
}
\source{
  Phelps, K. (1982).
  Use of the complementary log-log function to describe doseresponse
  relationships in insecticide evaluation field trials.
  \cr
  In R. Gilchrist (Ed.), \emph{Lecture Notes in Statistics, No. 14.
  GLIM.82: Proceedings of the International Conference on Generalized
  Linear Models}; Springer-Verlag.
}
\references{
  McCullagh P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models.}
  London: Chapman and Hall.

  Eva Cantoni and Elvezio Ronchetti (2001); JASA,  and \cr
  Eva Cantoni (2004); JSS, see \code{\link{glmrob}}
}

\examples{
data(carrots)
str(carrots)
plot(success/total ~ logdose, data = carrots, col = as.integer(block))
coplot(success/total ~ logdose | block, data = carrots)

## Classical glm
Cfit0 <- glm(cbind(success, total-success) ~ logdose + block,
             data=carrots, family=binomial)
summary(Cfit0)

## Robust Fit (see help(glmrob)) ....
}
\keyword{datasets}
