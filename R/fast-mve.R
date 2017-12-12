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

