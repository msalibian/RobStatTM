#' Robust R^2 coefficient of determination
#'
#' This function computes a robust version of the R^2 coefficient of determination.
#'
#' This function computes a robust version of the R^2 coefficient. It currently only
#' works for MM regression estimators computed with a rho function in Tukey's
#' bisquare family.
#'
#' @param RR2 the proportional difference in loss functions (a naive robust R^2 coefficient).
#' @param cc the tuning constant for the rho function
#'
#' @return An unbiased version of the robust R^2 coefficient of determination.
#'
#' @rdname INVTR2
#' @author Victor Yohai
#' @references \url{http://thebook}
#'
#' @export
INVTR2 <- function(RR2, family, cc) {

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
      hh <- function(v,family, cc, z) return(rho(v/z, family=family, cc, z)*dnorm(v) )
      ee <- 2*(integrate(hh, 0, cc[3]*z, family=family, cc=cc, z=z)$value+1-pnorm(cc[3]*z))
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










