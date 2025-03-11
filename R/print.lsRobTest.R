#' Print an lsRobTest Object
#'
#' @param x lmrobdetMM fitted model object
#' @param digits significant digits printed, default digits = 4
#' @param ... pass through parameters
#'
#' @returns print selected components of lmrobdetMM object
#' @export
#'
print.lsRobTest <- function(x, digits = 4, ...)
{
  cat("Test for least squares bias\n")
  if(x$test == "T1")
    cat("H0: normal regression error distribution\n")
  if(x$test == "T2")
    cat("H0: composite normal/non-normal regression error distribution\n")
  
  cat("\n")
  cat("Individual coefficient tests:\n")
  print(format(as.data.frame(x$coefs), digits = digits, ...))
  cat("\n")
  cat("Joint test for bias:\n")
  cat("Test statistic: ")
  cat(format(x$full$stat, digits = digits, ...))
  cat(" on ")
  cat(format(x$full$df, digits = digits, ...))
  cat(" DF, p-value: ")
  cat(format(x$full$p.value, digits = digits, ...))
  cat("\n")
  
  invisible(x)
}