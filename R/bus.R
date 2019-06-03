#' Bus data
#'
#' This data set corresponds to a study in automatic vehicle recognition.
#' Each of the 218 rows corresponds to a view of a bus silhouette, and contains
#' 18 attributes of the image. It was decided to exclude variable 9 and divide the
#' remaining variables by their MADN’s.
#'
#' Description: The following features were extracted from the silhouettes.
#' 1. compactness
#' 2. circularity
#' 3. distance circularity
#' 4. radius ratio
#' 5. principal axis aspect ratio
#' 6. maximum length aspect ratio
#' 7. scatter ratio
#' 8. elongatedness
#' 9. principal axis rectangularity
#' 10. maximum length rectangularity
#' 11. scaled variance along major axis
#' 12. scaled variance along minor axis
#' 13. scaled radius of gyration
#' 14. skewness about major axis
#' 15. skewness about minor axis
#' 16. kurtosis about minor axis
#' 17. kurtosis about major axis
#' 18. hollows ratio
#'
#' Format: Numeric, 218 rows and 18 columns.
#'
#' @docType data
#'
#' @usage data(bus)
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @source Hettich, S. and Bay, S.D. (1999), The UCI KDD Archive http://kdd.ics.uci.edu.
#' Irvine, CA: University of California, Department of Information and Computer Science.
#'
#' @examples
#' data(bus)
"bus"
