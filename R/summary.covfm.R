#' @export
summary.covfm <- function(object, ...) {
  chkDots(...)
  
  object.summary <- sapply(object, summary)
  oldClass(object.summary) <- "summary.covfm"
  object.summary
}

print.summary.covfm <- function(x, digits = max(3, getOption("digits") - 3),
                                     print.distance = FALSE, ...) {
  n.models  <- length(x)
  mod.names <- format(paste(names(x), ":", sep = ""), justify = "right")

  ## Need to add function calls for this to work
  # cat("\nCalls: \n")
  # for(i in 1:n.models) {
  #   cat(mod.names[i], ": ")
  #   print(x[[i]]$call)
  # }

  p <- dim(x[[1]]$cov)[1]
  i1 <- rep(seq(p), times = p)
  i2 <- rep(seq(p), each = p)

  cov.index <- paste("[", paste(i1, i2, sep = ","), "]", sep = "")
  cov.index <- matrix(cov.index, p, p)
  cov.index <- cov.index[row(cov.index) >= col(cov.index)]

  cov.unique <- t(sapply(x, function(u) u$cov[row(u$cov) >= col(u$cov)]))
  dimnames(cov.unique) <- list(mod.names, cov.index)

  cat("\nComparison of Covariance Estimates:\n")
  print(cov.unique, digits = digits, ...)
  
  center <- t(sapply(x, function(u) u$center))
  center.names <- names(x[[1]]$center)
  dimnames(center) <- list(mod.names, center.names)

  cat("\nComparison of Location Estimates:\n")
  print(center, digits = digits, ...)
  
  #if (all(sapply(x, function(u) u$corr) == FALSE)) {
  cov.index <- paste("[", paste(1:p, 1:p, sep = ","), "]", sep = "")
  
  cov.sd <- t(sapply(x, function(u) sqrt(diag(u$cov))))
  dimnames(cov.sd) <- list(mod.names, cov.index)
    
  cat("\nComparison of Scale Estimates:\n")
  print(cov.sd, digits = digits, ...)

  evals <- t(sapply(x, function(u) u$evals))
  eval.names <- names(x[[1]]$evals)
  dimnames(evals) <- list(mod.names, eval.names)

  cat("\nComparison of Eigenvalues: \n")
  print(evals, digits = digits, ...)

  have.dist <- sapply(x, function(u) !is.null(u$dist))
  if(print.distance && all(have.dist)) {
    dists <- t(sapply(x, function(u) u$dist))
    dimnames(dists) <- list(mod.names, names(x[[1]]$dist))
    cat("\nComparison of Mahalanobis Distances: \n")
    print(dists, digits = digits, ...)
  }

  invisible(x)
}