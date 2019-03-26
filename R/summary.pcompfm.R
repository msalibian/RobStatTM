#' @export
summary.pcompfm <- function(object, ...) {
  chkDots(...)
  
  object.summary <- sapply(object, summary)
  oldClass(object.summary) <- "summary.pcompfm"
  object.summary
}

#' @export
print.summary.pcompfm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  mod.names <- format(paste(names(x), ":", sep = ""), justify = "right")
  
  cat("Importance of Components:\n")
  
  cat("\nStandard deviation:\n")
  st.dev <- t(sapply(x, function(u) u$importance[1, ]))
  rownames(st.dev) <- mod.names
  print(st.dev, digits = digits, ...)
  
  cat("\nProportion of Variance:\n")
  prop.var <- t(sapply(x, function(u) u$importance[2, ]))
  rownames(prop.var) <- mod.names
  print(prop.var, digits = digits, ...)
  
  cat("\nCumulative Proportion\n")
  cumulative.prop <- t(sapply(x, function(u) u$importance[3, ]))
  rownames(cumulative.prop) <- mod.names
  print(cumulative.prop, digits = digits, ...)

  invisible(x)
}