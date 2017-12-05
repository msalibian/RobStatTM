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
INVTR2 <- function(RR2, cc) {

    TR2 <- function(R2, cc) {
    a <- Erhobic(cc)
    b <- Erhobic( cc*sqrt(1-R2))
    return( (b-a)/ (b*(1-a)) )
  }
  
  # compute E(rho(u,cc)), rho is the bisquare function
  Erhobic <- function(cc) {
    a0 <- 2*pnorm(cc)-1
    a2 <- (-2)*cc*dnorm(cc)+a0
    a4 <- (-2)*(cc^3)*dnorm(cc)+3*a2
    a6 <- (-2)*(cc^5)*dnorm(cc)+5*a4
    dd <- (a6/cc^6)+(3*a2/cc^2)-(3*a4/cc^4)+1-a0
    return( dd )
  }
  

  
  ff <- function(x, y, cc) return( TR2(x,cc) - y )
  
  
  aa <- TR2(.99999, cc)
  bb <- TR2(.00001, cc)
  if( RR2 > .99 ) R2 <- 1
  if ( RR2 < bb ) R2 <- 0 
  if( (RR2 <= .99) & (RR2 >= bb) ) 
    R2 <- uniroot(ff, c(bb/2, aa+((1-aa)/2)), y=RR2, cc=cc)$root #R2 <- uniroot(ff,c( .000000001,.999999999),y=RR2,cc=cc)$root
  return(R2)
}


# INVTR2(0.7367034,3.44)










