#' Tuning parameters for lmrob2
#'
#' This function sets tuning parameters for the MM-based Distance Constrained
#' Maximum Likelihood regression estimators computed by \code{lmrob2}.
#' 
#' There are 2 sets of tuning parameters: those related to the MM-estimator,
#' and those controlling the initial Pen~a-Yohai estimator. 
#'
#' @param seed \code{NULL}
#' @param tuning.chi tuning constant for the function used to compute the M-scale
#' for the S-estimator. For the estimator to be consistent it needs to 
#' be matched with the value of \code{bb} below. It defaults to 1.5477, which 
#' together with \code{bb = 0.5} yields an estimator with maximum breakdown point.
#' @param bb tuning constant (between 0 and 1) for the M-scale used to compute the initial S-estimator. It
#' determines the robusness (breakdown point) of the resulting MM-estimator, which is
#' \code{min(bb, 1-bb)}. Defaults to 0.5
#' @param tuning.psi tuning constant for the re-descending M-estimator. Its default
#' value (3.4434) returns an estimator with an asymptoti efficiency of 85% when errors
#' have a normal distribution
#' @param max.it maximum number of IRWLS iterations for the MM-estimator
#' @param refine.tol relative covergence tolerance for the S-estimator
#' @param rel.tol relative covergence tolerance for the IRWLS iterations for the MM-estimator
#' @param refine.PY number of refinement steps for the Pen~a-Yohai candidates
#' @param solve.tol relative tolerance for inversion
#' @param trace.lev positive values (increasingly) provide details on the progress of the MM-algorithm 
#' @param mts maximum number of samples
#' @param compute.rd logical value indicating whether robust leverage distances need to be computed. 
#' @param psi string specifying the type of loss function to be used.
#' @param corr.b logical value indicating whether a finite-sample correction should be applied 
#' to the M-scale parameter \code{bb}
#' @param split.type determines how categorical and continuous variables are split. See 
#' \link{\code{splitFrame}} in package \link{\code{robustbase}}. 
#' @param initial string specifying the initial value for the M-step of the MM-estimator. Valid
#' options are \code{'S'}, for an S-estimator and \code{'MS'} for an M-S estimator which is 
#' appropriate when there are categorical explanatory variables in the model.
#' @param prosac For \code{pyinit}, proportion of observations to remove based on PSCs. See \link{\code{pyinit}}.
#' @param clean.method For \code{pyinit}, how to clean the data based on large residuals. If 
#' \code{"threshold"}, all observations with scaled residuals larger than \code{C.res} will 
#' be removed, if \code{"proportion"}, observations with the largest \code{prop} residuals will 
#' be removed. See \link{\code{pyinit}}.
#' @param C.res See parameter \code{clean.method} above. See \link{\code{pyinit}}.
#' @param prop See parameter \code{clean.method} above. See \link{\code{pyinit}}.
#' @param py.nit Maximum number of iterations. See \link{\code{pyinit}}.
#' @param en.tol Relative tolerance for convergence.  See \link{\code{pyinit}}.
#' @param mscale.maxit Maximum number of iterations for the M-scale algorithm. See \link{\code{pyinit}}. 
#' @param mscale.tol Convergence tolerance for the M-scale algorithm. See \link{\code{pyinit}}.
#' @param mscale.rho.fun String indicating the loss function used for the M-scale. See \link{\code{pyinit}}.
#'
#' @return A list with the necessary tuning parameters. 
#' 
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' 
#' @seealso \code{\link{pyinit}}
#'
#' @examples
#' data(coleman)
#' m2 <- lmrob2(Y ~ ., data=coleman, control=lmrob2.control(refine.PY=50))
#'
#' @rdname lmrob2.control
#' @export
lmrob2.control <-  function(seed = NULL, tuning.chi = 1.5477, bb = 0.5, # 50% Breakdown point
                            tuning.psi = 3.4434, # 85% efficiency
                            max.it = 100, refine.tol = 1e-7, rel.tol = 1e-7,
                            refine.PY = 10, # no. of steps to refine PY candidates
                            solve.tol = 1e-7, trace.lev = 0, mts = 1000,
                            compute.rd = FALSE, psi = 'bisquare',
                            corr.b = TRUE, # for MMPY and SMPY
                            split.type = "f", # help(splitFrame, package='robustbase')
                            initial='S', #'S' or 'MS'
                            prosac = 0.5, clean.method = 'threshold', 
                            C.res = 2, prop = .2, py.nit = 20, en.tol = 1e-5, 
                            mscale.maxit = 50, mscale.tol = 1e-06, 
                            mscale.rho.fun = 'bisquare') {
  return(list(seed = as.integer(seed), psi=psi,
         tuning.chi=tuning.chi, bb=bb, tuning.psi=tuning.psi,
         max.it=max.it,
         refine.tol=refine.tol,
         corr.b = corr.b, refine.PY = refine.PY, 
         rel.tol=rel.tol,
         solve.tol=solve.tol, trace.lev=trace.lev, mts=mts,
         compute.rd=compute.rd, 
         split.type=split.type, 
         initial=initial, # method=method, subsampling=subsampling,
         prosac=prosac, clean.method=clean.method, C.res=C.res,
         prop=prop, py.nit=py.nit, en.tol=en.tol, mscale.maxit=mscale.maxit,
         mscale.tol=mscale.tol, mscale.rho.fun='bisquare'))
}

