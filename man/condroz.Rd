\name{condroz}
\alias{condroz}
\encoding{utf8}
\docType{data}
\title{ Condroz Data }
\description{
  Dataset with pH-value and Calcium content in soil samples, collected in different
  communities of the Condroz region in Belgium. The data pertain to a
  subset of 428 samples with a pH-value between 7.0 and 7.5.
}
\usage{data(condroz)}
\format{
  A data frame with 428 observations on the following 2 variables.
  \describe{
    \item{\code{Ca}}{Calcium content of the soil sample}
    \item{\code{pH}}{pH value of the soil sample}
  }
}
\details{
  For more information on the dataset, cf. Goegebeur et al. (2005).
}
\source{
  Hubert and Vandervieren (2006), p. 10.
  This dataset is also studied in Vandewalle et al. (2004).
}
\references{
  See also those for \code{\link{adjbox}}.

  Goegebeur, Y., Planchon, V., Beirlant, J., Oger, R. (2005). Quality
  Assesment of Pedochemical Data Using Extreme Value Methodology,
  Journal of Applied Science, 5, p. 1092-1102.

  Vandewalle, B., Beirlant, J., Hubert, M. (2004). A robust estimator of
  the tail index based on an exponential regression model, in Hubert,
  M., Pison G., Struyf, A. and S. Van Aelst, ed., Theory and
  Applications of Recent Robust Methods, Birkhäuser, Basel, p. 367-376.
}
\examples{
  adjbox(condroz$Ca)
}
\keyword{datasets}
