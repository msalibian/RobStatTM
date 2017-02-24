\name{adjbox}
\title{Plot an Adjusted Boxplot for Skew Distributions}
\alias{adjbox}
\alias{adjbox.default}
\alias{adjbox.formula}
\description{
  Produces boxplots adjusted for skewed distributions as proposed in
  Hubert and Vandervieren (2004).
}
\usage{
adjbox(x, \dots)

\method{adjbox}{formula}(formula, data = NULL, \dots, subset, na.action = NULL)

\method{adjbox}{default}(x, \dots, range = 1.5, doReflect = FALSE,
        width = NULL, varwidth = FALSE,
        notch = FALSE, outline = TRUE, names, plot = TRUE,
        border = par("fg"), col = NULL, log = "",
        pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
        horizontal = FALSE, add = FALSE, at = NULL)
}
\arguments{
    \item{formula}{a formula, such as \code{y ~ grp}, where \code{y} is a
    numeric vector of data values to be split into groups according to
    the grouping variable \code{grp} (usually a factor).}
  \item{data}{a data.frame (or list) from which the variables in
    \code{formula} should be taken.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used for plotting.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is to ignore missing
    values in either the response or the group.}
  \item{x}{for specifying data from which the boxplots are to be
    produced. Either a numeric vector, or a single list containing such
    vectors. Additional unnamed arguments specify further data
    as separate vectors (each corresponding to a component boxplot).
    \code{\link{NA}}s are allowed in the data.}
  \item{\dots}{For the \code{formula} method, named arguments to be passed to
    the default method.

    For the default method, unnamed arguments are additional data
    vectors (unless \code{x} is a list when they are ignored),
    and named arguments are arguments and graphical parameters to be
    passed to \code{\link{bxp}} in addition to the ones
    given by argument \code{pars} (and override those in \code{pars}).
  }
  \item{range}{this determines how far the plot whiskers extend out
    from the box, and is simply passed as argument \code{coef} to
    \code{\link{adjboxStats}()}.  If \code{range} is positive, the
    whiskers extend to the most extreme data point which is no more than
    \code{range} times the interquartile range from the box.  A value
    of zero causes the whiskers to extend to the data extremes.}
  \item{doReflect}{logical indicating if the MC should also be
    computed on the \emph{reflected} sample \code{-x}, and be averaged,
    see \code{\link{mc}}.}
  \item{width}{a vector giving the relative widths of the boxes making
    up the plot.}
  \item{varwidth}{if \code{varwidth} is \code{TRUE}, the boxes are
    drawn with widths proportional to the square-roots of the number
    of observations in the groups.}
  \item{notch}{if \code{notch} is \code{TRUE}, a notch is drawn in
    each side of the boxes.  If the notches of two plots do not
    overlap this is \sQuote{strong evidence} that the two medians differ
    (Chambers \emph{et al.}, 1983, p. 62).  See \code{\link{boxplot.stats}}
    for the calculations used.}
  \item{outline}{if \code{outline} is not true, the outliers are
    not drawn (as points whereas S+ uses lines).}% the argument name is most ugly but S+ compatible
  \item{names}{group labels which will be printed under each boxplot.}
  \item{boxwex}{a scale factor to be applied to all boxes.  When there
    are only a few groups, the appearance of the plot can be improved
    by making the boxes narrower.}
  \item{staplewex}{staple line width expansion, proportional to box
    width.}
  \item{outwex}{outlier line width expansion, proportional to box
    width.}
  \item{plot}{if \code{TRUE} (the default) then a boxplot is
    produced.  If not, the summaries which the boxplots are based on
    are returned.}
  \item{border}{an optional vector of colors for the outlines of the
    boxplots.  The values in \code{border} are recycled if the
    length of \code{border} is less than the number of plots.}
  \item{col}{if \code{col} is non-null it is assumed to contain colors
    to be used to colour the bodies of the box plots. By default they
    are in the background colour.}
  \item{log}{character indicating if x or y or both coordinates should
    be plotted in log scale.}
  \item{pars}{a list of (potentially many) more graphical parameters,
    e.g., \code{boxwex} or \code{outpch}; these are passed to
    \code{\link{bxp}} (if \code{plot} is true); for details, see there.}
  \item{horizontal}{logical indicating if the boxplots should be
    horizontal; default \code{FALSE} means vertical boxes.}
  \item{add}{logical, if true \emph{add} boxplot to current plot.}
  \item{at}{numeric vector giving the locations where the boxplots should
    be drawn, particularly when \code{add = TRUE};
    defaults to \code{1:n} where \code{n} is the number of boxes.}
}
\details{
  The generic function \code{adjbox} currently has a default method
  (\code{adjbox.default}) and a formula interface (\code{adjbox.formula}).

  If multiple groups are supplied either as multiple arguments or via a
  formula, parallel boxplots will be plotted, in the order of the
  arguments or the order of the levels of the factor (see
  \code{\link{factor}}).

  Missing values are ignored when forming boxplots.

  Extremes of the upper and whiskers of the adjusted boxplots are
  computed using the medcouple (\code{\link{mc}()}), a robust measure of
  skewness. For details, cf. TODO %% << FIXME
}
\value{
  A \code{\link{list}} with the following components:

  \item{stats}{a matrix, each column contains the extreme of the lower
    whisker, the lower hinge, the median, the upper hinge and the extreme of
    the upper whisker for one group/plot.  If all the inputs have the same
    class attribute, so will this component.}
  \item{n}{a vector with the number of observations in each group.}
  \item{coef}{a matrix where each column contains the lower and upper
    extremes of the notch.}
  \item{out}{the values of any data points which lie beyond the extremes
    of the whiskers.}
  \item{group}{a vector of the same length as out whose elements
    indicate to which group the outlier belongs.}
  \item{names}{a vector of names for the groups.}
}
\references{
%% Hubert, M. and Vandervieren, E. (2006)
%% \emph{An Adjusted Boxplot for Skewed Distributions},
%% Technical Report TR-06-11, KU Leuven, Section of Statistics, Leuven.
%% \url{http://wis.kuleuven.be/stat/robust/Papers/TR0611.pdf}

  Hubert, M. and Vandervieren, E. (2008).
  An adjusted boxplot for skewed distributions,
  \emph{Computational Statistics and Data Analysis} \bold{52}, 5186--5201.
}
\author{ R Core Development Team, slightly adapted by Tobias Verbeke }
\note{ The code and documentation only slightly modifies the code of
  \code{\link{boxplot.default}}, \code{boxplot.formula} and
  \code{\link{boxplot.stats}}
}
\seealso{The medcouple, \code{\link{mc}}; \code{\link{boxplot}}.
}
\examples{
if(require("boot")) {
 ### Hubert and Vandervieren (2006), p. 10, Fig. 4.
 data(coal, package = "boot")
 coaldiff <- diff(coal$date)
 op <- par(mfrow = c(1,2))
 boxplot(coaldiff, main = "Original Boxplot")
 adjbox(coaldiff, main  = "Adjusted Boxplot")
 par(op)
}

### Hubert and Vandervieren (2006), p. 11, Fig. 6. -- enhanced
op <- par(mfrow = c(2,2), mar = c(1,3,3,1), oma = c(0,0,3,0))
with(condroz, {
 boxplot(Ca, main = "Original Boxplot")
 adjbox (Ca, main = "Adjusted Boxplot")
 boxplot(Ca, main = "Original Boxplot [log]", log = "y")
 adjbox (Ca, main = "Adjusted Boxplot [log]", log = "y")
})
mtext("'Ca' from data(condroz)",
      outer=TRUE, font = par("font.main"), cex = 2)
par(op)
}
\keyword{hplot}
