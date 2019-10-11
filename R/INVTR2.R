#' Robust R^2 coefficient of determination
#'
#' This function computes a robust version of the R^2 coefficient of determination.
#' It is used internally by \code{\link{lmrobdetMM}},
#' and not meant to be used directly.
#'
#' This function computes a robust version of the R^2 coefficient.
#' It is used internally by \code{\link{lmrobdetMM}},
#' and not meant to be used directly.
#'
#' @param RR2 the proportional difference in loss functions (a naive robust R^2 coefficient).
#' @param family family string specifying the name of the family of loss function to be used (current valid
#' options are "bisquare", "opt" and "mopt").
#' @param cc tuning parameters to be computed according to efficiency and / or breakdown
#' considerations. See \link{lmrobdet.control}, \link{bisquare}, \link{mopt}
#' and \link{opt}.

#'
#' @return An unbiased version of the robust R^2 coefficient of determination.
#'
#' @rdname INVTR2
#' @author Victor Yohai, \email{victoryohai@gmail.com}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @export
INVTR2 <- function(RR2, family, cc) {

  hh <- function(v, family, cc, z) return( rho(v/z, family=family, cc)*dnorm(v) )

  TR2 <- function(R2, family, cc) {
    a <- Erhobic(family, cc, 1)
    b <- Erhobic(family, cc, sqrt(1-R2))
    return( (b-a)/ (b*(1-a)) )
  }

  # compute E(rho(u,cc)), rho is the bisquare function
  Erhobic <- function(family, cc, zz) {
    if( family == 'bisquare') {
      dd <- cc * zz
      a0 <- 2*pnorm(dd)-1
      a2 <- (-2)*dd*dnorm(dd)+a0
      a4 <- (-2)*(dd^3)*dnorm(dd)+3*a2
      a6 <- (-2)*(dd^5)*dnorm(dd)+5*a4
      ee <- (a6/dd^6)+(3*a2/dd^2)-(3*a4/dd^4)+1-a0
    } else {
      ee <- 2*(integrate(hh, 0, cc[3]*zz, family=family, cc=cc, z=zz)$value+1-pnorm(cc[3]*zz))
    }
    return( ee )
  }

  ff <- function(x, y, family, cc) return( TR2(x, family, cc) - y )
  aa <- TR2(.99999, family, cc)
  bb <- TR2(.00001, family, cc)
  if( RR2 > .99 ) R2 <- 1
  if ( RR2 < bb ) R2 <- 0
  if( (RR2 <= .99) & (RR2 >= bb) )
    R2 <- uniroot(ff, c(bb/2, aa+((1-aa)/2)), y=RR2, family=family, cc=cc)$root #R2 <- uniroot(ff,c( .000000001,.999999999),y=RR2,cc=cc)$root
  return(R2)
}


# INVTR2(0.7367034,3.44)










