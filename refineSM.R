refine.sm <- function(x, y, initial.beta, initial.scale, k=50, 
                      conv=1, b, cc, step='M') {
  # does "k" IRWLS refining steps from "initial.beta"
  #
  # if "initial.scale" is present, it's used, o/w the MAD is used
  # k = number of refining steps
  # conv = 0 means "do k steps and don't check for convergence"
  # conv = 1 means "stop when convergence is detected, or the
  #                 maximum number of iterations is achieved"
  # b and cc = tuning constants of the equation
  # step = 'M' means M-IRWLS iterations (scale is not updated)
  # step = 'S' means S-IRWLS iterations (scale is updated)
  # 
  
  
  n <- dim(x)[1]
  # p <- dim(x)[2]
  
  res <- as.vector( y - x %*% initial.beta )
  
  if( missing( initial.scale ) ) {
    initial.scale <- scale <- median(abs(res))/.6745
  } else {
    scale <- initial.scale
  }
  
  beta <- initial.beta
  
  
  converged <- FALSE
  
  # lower.bound <- median(abs(res))/cc
  
  
  for(i in 1:k) {
    # do one step of the iterations to solve for the scale
    scale.super.old <- scale
    #lower.bound <- median(abs(res))/1.56
    if(step=='S') {
      scale <- sqrt( scale^2 * mean( rho( res / scale, cc ) ) / b     )
      # scale <- mscale(res, tol=1e-7, delta=b, max.it=500, tuning.chi=cc)
    }
    # now do one step of IRWLS with the "improved scale"
    weights <- f.w( res/scale, cc )
    # W <- matrix(weights, n, p)
    xw <- x * sqrt(weights) # sqrt(W)
    yw <- y *   sqrt(weights)
    beta.1 <- our.solve( t(xw) %*% xw ,t(xw) %*% yw )
    if(any(is.na(beta.1))) { beta.1 <- initial.beta
    scale <- initial.scale
    break
    }
    if( (conv==1) )
    {
      # check for convergence
      if( norm.sm( beta - beta.1 ) / norm.sm(beta) < 1e-7 ) { # magic number alert!!!
        converged <- TRUE
        break
      }
    }
    res <- as.vector( y - x %*% beta.1 )
    beta <- beta.1
    # print(as.vector(t(x) %*% rhoprime(res/scale, cc))/n)
    # print(scale)
  }
  
  # res <- as.vector( y - x %*% beta )
  # get the residuals from the last beta
  return(list(beta.rw = beta.1, scale.rw = scale, converged=converged))
  
}

norm.sm <- function(x) sqrt(sum(x^2))



our.solve <- function(a,b) {
  a <- qr(a)
  da <- dim(a$qr)
  if(a$rank < (p <- da[2]))
    return(NA)
  else qr.coef(a, b)
}


## Weight function
f.w <- function(u, cc) {
  tmp <- (1 - (u/cc)^2)^2
  tmp[abs(u/cc) > 1] <- 0
  return(tmp)
}


