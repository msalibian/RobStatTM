#' Minimum Volume Ellipsoid covariance estimator
#'
#' This function uses a fast algorithm to compute the Minimum Volume
#' Ellipsoid (MVE) for multivariate location and scatter.
#'
#' This function computes the Minimum Volume
#' Ellipsoid (MVE) for multivariate location and scatter, using a
#' fast algorithm related to the fast algorithm for S-regression
#' estimators (see \code{\link[robustbase]{lmrob}}).
#'
#' @param x data matrix (n x p) with cases stored in rows.
#' @param nsamp number of random starts for the iterative algorithm, these
#' are constructed using subsamples of the data.
#'
#' @return A list with the following components:
#' \item{center}{a vector with the robust multivariate location estimate}
#' \item{cov}{a matrix with the robust covariance / scatter matrix estimate}
#' \item{scale}{A scalar that equals the median of the mahalanobis distances of the
#' data to the \code{center}, multiplied by the determinant of the covariance matrix
#' to the power 1/p}
#' \item{best}{Indices of the observations that correspond to the MVE estimator}
#' \item{nsamp}{Number of random starts used for the iterative algorithm}
#' \item{nsing}{Number of random subsamples (among the \code{nsamp} attempted)
#' that failed (resulting in singular initial values)}
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(bus)
#' X0 <- as.matrix(bus)
#' X1 <- X0[,-9]
#' tmp <- fastmve(X1)
#' round(tmp$cov[1:10, 1:10], 3)
#' tmp$center
#'
#' @export
fastmve <- function(x, nsamp=500) {
  n <- nrow(x)
  p <- ncol(x)
  n2 <- floor(n/2)
  nind <- p +1
  tmp <- .C('r_fast_mve', as.double(x),
            as.integer(n), as.integer(p), as.integer(nsamp),
            nsing = as.integer(0), ctr = as.double(rep(0,p)),
            cov = as.double(rep(0,p*p)),
            scale = as.double(0), best=as.integer(rep(0,n)),
            as.integer(nind), as.integer(n2), PACKAGE='RobStatTM')
  mve.cov <- matrix(tmp$cov, p, p)
  return(list(center= tmp$ctr, cov=mve.cov, scale=tmp$scale,
              best=tmp$best[1:floor(n/2)],
              nsamp=nsamp, nsing = tmp$nsing))
}

