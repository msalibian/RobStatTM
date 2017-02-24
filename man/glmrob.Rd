\name{glmrob}
\alias{glmrob}
\title{Robust Fitting of Generalized Linear Models}
\encoding{utf8}
\description{
  \code{glmrob} is used to fit generalized linear models by robust
  methods.  The models are specified by giving a symbolic description of
  the linear predictor and a description of the error distribution.
  Currently, robust methods are implemented for \code{\link{family} =
  binomial},  \code{ = poisson}, \code{ = Gamma} and \code{ = gaussian}.
}

\usage{
glmrob(formula, family, data, weights, subset, na.action,
       start = NULL, offset, method = c("Mqle", "BY", "WBY", "MT"),
       weights.on.x = c("none", "hat", "robCov", "covMcd"), control = NULL,
       model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, trace.lev = 0, ...)
}
\arguments{
  \item{formula}{a \code{\link{formula}}, i.e., a symbolic description
    of the model to be fit (cf. \code{\link{glm}} or \code{\link{lm}}).}

  \item{family}{a description of the error distribution and link function to
    be used in the model.  This can be a character string naming a
    family function, a family \code{\link{function}} or the result of a
    call to a family function.  (See \code{\link{family}} for details of
    family functions.)}

  \item{data}{an optional data frame containing the variables
    in the model.  If not found in \code{data}, the variables are taken
    from \code{environment(formula)}, typically the environment from
    which \code{glmrob} is called.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process.}

  \item{subset}{an optional vector specifying a subset of observations to be
    used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting in \code{\link{options}}. The
    \dQuote{factory-fresh} default is \code{\link{na.omit}}.}

  \item{start}{starting values for the parameters in the linear
    predictor.  Note that specifying \code{start} has somewhat different
    meaning for the different \code{method}s.  Notably, for \code{"MT"},
    this skips the expensive computation of initial estimates via sub
    samples, but needs to be \emph{robust} itself.}

  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.}

  \item{method}{a character string specifying the robust fitting
    method.  The details of method specification are given below.}

  \item{weights.on.x}{
    a character string (can be abbreviated), a \code{\link{function}} or
    \code{\link{list}} (see below), or a numeric vector of length
    \code{n}, specifying how points (potential outliers) in x-space are
    downweighted.  If \code{"hat"},
    weights on the design of the form \eqn{\sqrt{1-h_{ii}}} are used,
    where \eqn{h_{ii}} are the diagonal elements of the hat matrix.  If
    \code{"robCov"}, weights based on the robust Mahalanobis distance of the
    design matrix (intercept excluded) are used where the covariance
    matrix and the centre is estimated by \code{\link[MASS]{cov.rob}}
    from the package \pkg{MASS}.\cr
    Similarly,  if \code{"covMcd"}, robust weights are computed using
    \code{\link{covMcd}}.  The default is \code{"none"}.

    If \code{weights.on.x} is a \code{\link{function}}, it is called
    with arguments \code{(X, intercept)} and must return an n-vector of
    non-negative weights.

    If it is a \code{\link{list}}, it must be of length one, and as
    element contain a function much like \code{\link{covMcd}()} or
    \code{\link[MASS]{cov.rob}()} (package \pkg{MASS}), which computes
    multivariate location and \dQuote{scatter} of a data matrix \code{X}.
  }
  \item{control}{a list of parameters for controlling the fitting process.
    See the documentation for \code{\link{glmrobMqle.control}} for
    details.}

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{x, y}{logical values indicating whether the response
    vector and model matrix used in the fitting process should be
    returned as components of the returned value.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}
  \item{trace.lev}{logical (or integer) indicating if intermediate results
    should be printed; defaults to \code{0} (the same as \code{FALSE}).}
  \item{\dots}{arguments passed to \code{\link{glmrobMqle.control}} when
    \code{control} is \code{NULL} (as per default).}
}

\details{
  \code{method="model.frame"} returns the \code{\link{model.frame}()},
  the same as \code{\link{glm}()}.
  \cr
  \code{method="Mqle"} fits a generalized linear
  model using Mallows or Huber type robust estimators, as described in
  Cantoni and Ronchetti (2001) and Cantoni and Ronchetti (2006).  In
  contrast to the implementation
  described in Cantoni (2004), the pure influence algorithm is
  implemented.
  \cr
  \code{method="WBY"} and \code{method="BY"},
  available for logistic regression (\code{family = binomial}) only, call
  \code{\link{BYlogreg}(*, initwml= . )} for the (weighted) Bianco-Yohai
  estimator, where \code{initwml} is true for \code{"WBY"}, and false
  for \code{"BY"}.
  \cr
  \code{method="MT"}, currently only implemented for \code{family = poisson},
  computes an \dQuote{[M]-Estimator based on [T]ransformation},
  by Valdora and Yohai (2013).

  \code{weights.on.x= "robCov"} makes sense if all explanatory variables
  are continuous.

  In the cases,where \code{weights.on.x} is \code{"covMcd"} or
  \code{"robCov"}, or list with a \dQuote{robCov} function, the
  mahalanobis distances \code{D^2} are computed with respect to the
  covariance (location and scatter) estimate, and the weights are
  \code{1/sqrt(1+ pmax.int(0, 8*(D2 - p)/sqrt(2*p)))},
  where \code{D2 = D^2} and \code{p = ncol(X)}.
}

\value{
  \code{glmrob} returns an object of class \code{"glmrob"} and is also
  inheriting from \code{\link{glm}}.
  \cr
  The \code{\link{summary}} method, see \code{\link{summary.glmrob}}, can
  be used to obtain or print a summary of the results.
  \cr
  The generic accessor functions \code{\link{coefficients}},
  \code{effects}, \code{fitted.values} and \code{residuals} (see
  \code{\link{residuals.glmrob}}) can be used to extract various useful
  features of the value returned by \code{glmrob()}.

  An object of class \code{"glmrob"} is a list with at least the
  following components:
  \item{coefficients}{a named vector of coefficients}
  \item{residuals}{the \emph{working} residuals, that is the (robustly
  \dQuote{huberized}) residuals in the final iteration of the IWLS fit.}
  \item{fitted.values}{the fitted mean values, obtained by transforming
    the linear predictors by the inverse of the link function.}
  \item{w.r}{robustness weights for each observations; i.e.,
    \code{residuals} \eqn{\times}{*} \code{w.r} equals the psi-function of the
    Preason's residuals.}
  \item{w.x}{weights used to down-weight observations based on the
    position of the observation in the design space.}
  \item{dispersion}{robust estimation of dispersion paramter if appropriate}
  \item{cov}{the estimated asymptotic covariance matrix of the estimated
    coefficients.}
  \item{tcc}{the tuning constant c in Huber's psi-function.}
  \item{family}{the \code{\link{family}} object used.}
  \item{linear.predictors}{the linear fit on link scale.}
  \item{deviance}{NULL; Exists because of compatipility reasons.}
  \item{iter}{the number of iterations used by the influence algorithm.}
  \item{converged}{logical. Was the IWLS algorithm judged to have converged?}
  \item{call}{the matched call.}
  \item{formula}{the formula supplied.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{data}{the \code{data argument}.}
  \item{offset}{the offset vector used.}
  \item{control}{the value of the \code{control} argument used.}
  \item{method}{the name of the robust fitter function used.}
  \item{contrasts}{(where relevant) the contrasts used.}
  \item{xlevels}{(where relevant) a record of the levels of the factors
    used in fitting.}

%% FIXME: This is for glm() -- but *not* (yet ??) for glmrob()
%% -----  should we change?
  % If a \code{\link{binomial}} \code{glm} model was specified by giving a
  % two-column response, the weights returned by \code{prior.weights} are
  % the total numbers of cases (multipied by the supplied case weights) and
  % the component \code{y} of the result is the proportion of successes.
}

\references{
  Eva Cantoni and Elvezio Ronchetti (2001)
  Robust Inference for Generalized Linear Models.
  \emph{JASA} \bold{96} (455), 1022--1030.

  Eva Cantoni (2004)
  Analysis of Robust Quasi-deviances for Generalized Linear Models.
  \emph{Journal of Statistical Software}, \bold{10},
  \url{http://www.jstatsoft.org/v10/i04}

  Eva Cantoni and Elvezio Ronchetti (2006)
  A robust approach for skewed and heavy-tailed outcomes in the analysis
  of health care expenditures.
  \emph{Journal of Health Economics} \bold{25}, 198--213.

  S. Heritier, E. Cantoni, S. Copt, M.-P. Victoria-Feser (2009)
  \emph{Robust Methods in Biostatistics}. Wiley Series in Probability
  and Statistics.

  Marina Valdora and Víctor J. Yohai (2013)
  Robust estimators for Generalized Linear Models.  In progress.
}

\author{Andreas Ruckstuhl ("Mqle") and Martin Maechler}
%%\note{ }

\seealso{
  \code{\link{predict.glmrob}} for prediction;
  \code{\link{glmrobMqle.control}}
}
\examples{
## Binomial response --------------
data(carrots)

Cfit1 <- glm(cbind(success, total-success) ~ logdose + block,
             data = carrots, family = binomial)
summary(Cfit1)

Rfit1 <- glmrob(cbind(success, total-success) ~ logdose + block,
                family = binomial, data = carrots, method= "Mqle",
                control= glmrobMqle.control(tcc=1.2))
summary(Rfit1)

Rfit2 <- glmrob(success/total ~ logdose + block, weights = total,
                family = binomial, data = carrots, method= "Mqle",
                control= glmrobMqle.control(tcc=1.2))
coef(Rfit2)  ## The same as Rfit1


## Binary response --------------
data(vaso)

Vfit1 <- glm(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso)
coef(Vfit1)

Vfit2 <- glmrob(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso,
                method="Mqle", control = glmrobMqle.control(tcc=3.5))
coef(Vfit2) # c = 3.5 ==> not much different from classical
## Note the problems with  tcc <= 3 %% FIXME algorithm ???

Vfit3 <- glmrob(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso,
                method= "BY")
coef(Vfit3)## note that results differ much.
## That's not unreasonable however, see Kuensch et al.(1989), p.465

## Poisson response --------------
data(epilepsy)

Efit1 <- glm(Ysum ~ Age10 + Base4*Trt, family=poisson, data=epilepsy)
summary(Efit1)

Efit2 <- glmrob(Ysum ~ Age10 + Base4*Trt, family = poisson,
                data = epilepsy, method= "Mqle",
                control = glmrobMqle.control(tcc= 1.2))
summary(Efit2)

## 'x' weighting:
(Efit3 <- glmrob(Ysum ~ Age10 + Base4*Trt, family = poisson,
 	         data = epilepsy, method= "Mqle", weights.on.x = "hat",
		 control = glmrobMqle.control(tcc= 1.2)))

try( # gives singular cov matrix: 'Trt' is binary factor -->
     # affine equivariance and subsampling are problematic
Efit4 <- glmrob(Ysum ~ Age10 + Base4*Trt, family = poisson,
                data = epilepsy, method= "Mqle", weights.on.x = "covMcd",
                control = glmrobMqle.control(tcc=1.2, maxit=100))
)

##--> See  example(possumDiv)  for another  Poisson-regression


### -------- Gamma family -- data from example(glm) ---

clotting <- data.frame(
            u = c(5,10,15,20,30,40,60,80,100),
         lot1 = c(118,58,42,35,27,25,21,19,18),
         lot2 = c(69,35,26,21,18,16,13,12,12))
summary(cl <- glm   (lot1 ~ log(u), data=clotting, family=Gamma))
summary(ro <- glmrob(lot1 ~ log(u), data=clotting, family=Gamma))

clotM5.high <- within(clotting, { lot1[5] <- 60 })
op <- par(mfrow=2:1, mgp = c(1.6, 0.8, 0), mar = c(3,3:1))
plot( lot1  ~ log(u), data=clotM5.high)
plot(1/lot1 ~ log(u), data=clotM5.high)
par(op)
## Obviously, there the first observation is an outlier with respect to both
## representations!

cl5.high <- glm   (lot1 ~ log(u), data=clotM5.high, family=Gamma)
ro5.high <- glmrob(lot1 ~ log(u), data=clotM5.high, family=Gamma)
with(ro5.high, cbind(w.x, w.r))## the 5th obs. is downweighted heavily!

plot(1/lot1 ~ log(u), data=clotM5.high)
abline(cl5.high, lty=2, col="red")
abline(ro5.high, lwd=2, col="blue") ## result is ok (but not "perfect")

%% FIXME: Need work -- option of *starting* from
%% -----  see Andreas' ~/R/MM/Pkg-ex/robustbase/glmrob-gamma-ARu.R

% ## a "regular outlier" in the middle :
% clotM4.3 <- within(clotting, { lot1[4] <- 1000 })
% ## .. not even this one works : ... need *robust* start ?!
% try(cl4.3 <- glm   (lot1 ~ log(u), data=clotM4.3, family=Gamma))
% try(ro4.3 <- glmrob(lot1 ~ log(u), data=clotM4.3, family=Gamma))

% ## The new option to start from "lmrobMM" --- not yet ok either
% try(
% ro4.3 <- glmrob(lot1 ~ log(u), data=clotM4.3, family=Gamma,
%                 start = "lmrobMM")
% )
% ## summary(ro4.3)

%% TODO the "same" with lot2 :
%% summary(glm(lot2 ~ log(u), data=clotting, family=Gamma))
}
\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
