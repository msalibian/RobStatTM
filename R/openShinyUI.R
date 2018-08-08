#' Open the Shiny interface for the package
#'
#' This function opens the Shiny interface for the package.
#'
#' @author Gregory Brownson \email{gsb25@uw.edu}
#' @references \url{http://thebook}
#' 
#' @export
openShinyUI <- function() {
  appDir <- system.file("shiny-ui", package = "RobStatTM")
  if (appDir == "") {
    stop("Could not find Shiny UI directory. Try re-installing `RobStatTM`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}  
