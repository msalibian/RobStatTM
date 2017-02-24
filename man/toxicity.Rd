\name{toxicity}
\encoding{utf8}
\alias{toxicity}
\docType{data}
\title{Toxicity of Carboxylic Acids Data}
\description{
  The aim of the experiment was to predict the toxicity of carboxylic acids on
  the basis of several molecular descriptors.
}
\usage{data(toxicity)}
\format{
  A data frame with 38 observations on the following 10 variables which are
  attributes for carboxylic acids:
  \describe{
    \item{\code{toxicity}}{aquatic toxicity, defined as
      \eqn{\log(IGC_{50}^{-1})}{log(IGC50^(-1))}; typically the \dQuote{response}.}
    \item{\code{logKow}}{\eqn{log Kow}, the partition coefficient}
    \item{\code{pKa}}{pKa: the dissociation constant}
    \item{\code{ELUMO}}{\bold{E}nergy of the \bold{l}owest
      \bold{u}noccupied \bold{m}olecular \bold{o}rbital}
    \item{\code{Ecarb}}{Electrotopological state of the \bold{carb}oxylic group}
    \item{\code{Emet}}{Electrotopological state of the \bold{met}hyl group}
    \item{\code{RM}}{Molar refractivity}
    \item{\code{IR}}{Refraction index}
    \item{\code{Ts}}{Surface tension}
    \item{\code{P}}{Polarizability}
  }
}
% \details{
% }
\source{
  The website accompanying the MMY-book:
  \url{http://www.wiley.com/legacy/wileychi/robust_statistics}
}
\references{
  Maguna, F.P., Núñez, M.B., Okulik, N.B. and Castro, E.A. (2003)
  Improved QSAR analysis of the toxicity of aliphatic carboxylic acids;
  \emph{Russian Journal of General Chemistry} \bold{73}, 1792--1798.
}
\examples{
data(toxicity)
summary(toxicity)
plot(toxicity)
plot(toxicity ~ pKa, data = toxicity)

## robustly scale the data (to scale 1) using Qn
(scQ.tox <- sapply(toxicity, Qn))
scTox <- scale(toxicity, center = FALSE, scale = scQ.tox)
csT <- covOGK(scTox, n.iter = 2,
              sigmamu = s_Qn, weight.fn = hard.rejection)
as.dist(round(cov2cor(csT$cov), 2))
}
\keyword{datasets}
