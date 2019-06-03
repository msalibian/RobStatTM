#' Algae data
#'
#' Each row of the data set is a set of 90 measurements at a river in some place in Europe.
#' There are 11 predictors. The response is the logarithm of the abundance of a certain
#' class of algae. Description: The columns are:
#' 1. season, categorical  (1,2,3,4 for winter, spring, summer and autumn)
#' 2. river size (categorical) (1,2,3 for small, medium and large)
#' 3. fluid velocity (categorical) (1,2,3 for low, medium and high)
#' 4-11 (numerci): content of nitrogen in the form of nitrates, nitrites and ammonia, and other
#' chemical compounds.
#' Col. 12 ia the response:  abundance of a type of algae (type 6 in the complete file). For
#' simplicity we deleted the rows with missing values and took the logarithm of the response.
#'
#' Format 90 rows, 12 columns (3 categorical, 9 numeric)
#'
#' @docType data
#'
#' @usage data(algae)
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @references References go here.
#'
#' @source Hettich, S. and Bay, S.D. (1999), The UCI KDD Archive http://kdd.ics.uci.edu.
#' Irvine, CA: University of California, Department of Information and Computer Science.
#'
#' @examples
#' data(algae)
"algae"
