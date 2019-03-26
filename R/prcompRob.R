#' Robust Principal Components Cont'd
#'
#' This function uses the pcaRobS function to compute all principal components while
#' behaving similarly to the prcomp function
#' 
#' @export prcompRob
#' @aliases prcompRob
#' @rdname prcompRob
#' 
#' @param x data matrix with observations in rows
#' @param rank. Maximal number of principal components to be used (optional)
#' @param delta.scale "delta" parametor of the scale M-estimator (default = 0.5)
#' @param max.iter maximum number of iterations (default = 100)
#' 
#' @return
#' \item{sdev}{the standard deviation of the principal components}
#' \item{rotation}{matrix containing the factor loadings}
#' \item{x}{matrix containing the rotated data}
#' \item{center}{the centering used}
#' 
#' @author Gregory Brownson, \email{gregory.brownson@gmail.com}
#'
#' @examples
#' data(wine)
#' 
#' p.wine <- prcompRob(wine)
#' summary(p.wine)
#' 
#' ## Choose only 5
#' p5.wine <- prcompRob(wine, rank. = 5)
#' summary(p5.wine)
#' 

prcompRob <- function(x, rank. = NULL, delta.scale = 0.5, max.iter = 100L) {
  ncomp <- if (is.null(rank.)) {
             ncol(x)
           } else {
             rank.
           }
  
  X = pcaRobS(x, ncomp = ncomp, desprop = 1.0, deltasca = delta.scale, maxit = max.iter)
  
  # X.scaled <- scale(data, center = center, scale = scale)
  
  n <- ncol(X$eigvec)
  
  pc.index <- sapply(1:n, function(i) paste("PC", i))
  
  z <- list()
  
  z$sdev     <- sqrt(diag(var(X$repre)))
  z$rotation <- X$eigvec
  colnames(z$rotation) <- pc.index
  
  z$center <- X$mu
  z$x      <- X$repre
  
  class(z) <- "prcompRob"
  
  z
}

#' @export
print.prcompRob <- function(x, print.x = FALSE, ...) {
    cat("Standard deviations:\n")
    print(x$sdev, ...)
    
    cat("\nEigenvectors:\n")
    print(x$rotation, ...)
    
    if (print.x && length(x$x)) {
        cat("\nRotated data:\n")
        print(x$x, ...)
    }
    
    invisible(x)
}

#' @export
summary.prcompRob <- function(object, ...) {
  chkDots(...)
  
  vars <- object$sdev^2
  vars <- vars/sum(vars)
  
  importance <- rbind("Standard deviation" = object$sdev,
                      "Proportion of Variance" = round(vars, 5),
                      "Cumulative Proportion" = round(cumsum(vars), 5))
  
  colnames(importance) <- colnames(object$rotation)
  
  object$importance <- importance
  
  class(object) <- "summary.prcompRob"
  object
}

#' @export
print.summary.prcompRob <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("Importance of components:\n")
    print(x$importance, digits = digits, ...)
    invisible(x)
}