#' @title Test for Least Squares Bias Using Robust MM Regressions
#'
#' @param object An MM regression fitted model whose class is *lmrobdetMM*.
#' @param test A character vector indicating which of two type of tests "T" or 
#' "T0: are used, with type "T" the default.
#' @param ... Pass through parameters
#'
#' @returns A list with component names coefs, full, test, efficiency
#' 
#' @details The original version of \code{lsRobTestMM} is the \code{lsRobTest}
#' in the package *robust*. The function \code{lsRobTest} had options *T1* and
#' *T2*. However, we only recommend *T2*, and deprecate *T1*. Accordingly we
#' use *T* for the former *T2*, and use *T0* for the former *T1*, and we
#' deprecate *T0*.
#' 
#' The *coefs* component of the list is a matrix whose row names are
#' the names of the regression predictor variables, and whose columns *LS*, 
#' *Robust*, *Delta*, *Std.error*, *t-stat*, *p-value* contain respectively,
#' the least squares and robust coefficient estimates, the differences in the
#' coefficient estimates, the standard errors of the differences, the resulting
#' t-statistic values, and the resulting z-test p-values.
#' 
#' The *full* component of the list is itself a list with components the full
#' model quadratic form chi-squared statistic value (*stat*), the degrees of 
#' freedom (*df*), and the full model p value (*p.value*).
#' 
#' The *test* component of the list is a character value indicating which of the 
#' tests *T* and *T0* has been computed.
#' 
#' The *efficiency* component of the list is *NULL* when test *T* has been used,
#' and is equal to the normal distribution efficiency of the *lmrobdetMM*
#' estimate when test *T0* has been used. 
#' 
#' @author Kjell Konis
#' 
#' @export
#'
#' @examples
#' args(lsRobTestMM)
lsRobTestMM <- function(object, test = c("T", "T0"), ...)
{
  # require(RobStatTM)
  test <- match.arg(test)
  
  # family: one of bisquare, opt and mopt
  family <- object$control$family
  tune <- object$control$tuning.psi
  
  # the prefix probably can be removed when added into RobStatTM
  eff = computeGaussianEfficiencyFromFamily(family, tune)
  
  if(is.null(object$weights)) {
    LS <- lm(formula(object$terms), data = object$model)
  } else {
    LS <- lm(formula(object$terms), data = object$model, weights = object$weights)
  }
  
  rmm <- residuals(object)
  rls <- residuals(LS)
  rob.sigma <- object$scale
  
  # require(robust) # probably not needed anymore, check later
  # tune <- lmRob.effvy(eff, ipsi) # control$tuning.psi
  # rw <- object$T.M.weights 
  
  rw <- object$rweights
  
  X <- model.matrix(object)
  n <- nrow(X)
  p <- ncol(X)
  
  if (!is.null(object$weights)) {
    X <- X * sqrt(object$weights)
  }
  
  V <- (t(rw * X) %*% X) / sum(rw) 
  V.inv <- solve(V)
  
  if(test == "T0") {
    d <- mean(rhoprime2(rmm/rob.sigma, family=family,cc=tune))
    tau <- n * mean( (rhoprime(rmm/rob.sigma, family=family,cc=tune)/d)^2 ) / (n - p)
    mat <- (1 - eff)/n * tau * rob.sigma^2 * V.inv 
  }
  
  if(test == "T") {
    d <- mean(rhoprime2(rmm/rob.sigma, family=family,cc=tune))
    delta2 <- mean( (rls - (rob.sigma * rhoprime(rmm/rob.sigma, family=family,cc=tune)) / d)^2 )
    mat <- delta2 / n * V.inv
  }
  
  brob <- coef(object)
  coef.names <- names(brob)
  bls <- coef(LS)
  x <- bls - brob
  
  if(attributes(object$terms)$intercept) {
    brob <- brob[-1]
    bls <- bls[-1]
    x <- x[-1]
    mat <- mat[-1, -1, drop = FALSE]
    coef.names <- coef.names[-1]
  }
  
  se <- sqrt(diag(mat))
  uniV <- x / se
  coefs <- cbind(bls, brob, x, se, uniV, 2*pnorm(-abs(uniV)))
  dimnames(coefs) <- list(coef.names, c("LS", "Robust", "Delta", "Std. Error", "Stat", "p-value"))
  T <- drop(t(x) %*% solve(mat) %*% x)
  
  ans <- list(coefs = coefs,
              full = list(stat = T, df = length(x), p.value = 1 - pchisq(T, length(x))),
              test = test,
              efficiency = eff)
  
  # class here is not changed yet, also the postfix in the name of print
  oldClass(ans) <- "lsRobTest"
  ans
}