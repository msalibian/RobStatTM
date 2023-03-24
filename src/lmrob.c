/* -*- mode: c; kept-new-versions: 40; kept-old-versions: 20 -*-
 * Indentation (etc) style:  C-c . gnu */

/* file lmrob.c
 * was	roblm/src/roblm.c - version 0.6	 by Matias Salibian-Barreras

 * Includes the fast-s algorithm
 */

/* Robust MM regression estimates *
 * ------------------------------ */

/* comment code *
 *
 * adapt other sampler <<<<<<<<<<  R's random number generator !!!!

 * replace abort for too many singular resamples by
 * returning the number of singular ones
 */

/* MM:
   -  Done:  fast_s[_large_n]() both had FIXED seed (= 37),
	     and effectively discarded the seed_rand argument below
   -  Done:  drop 'register' : today's compilers do optimize well!
   -  Done:  using Calloc() / Free() instead of malloc()/free()
*/

/* kollerma:
   Added alternative psi functions callable via psifun, chifun and
   wgtfun. ipsi is used to distinguish between the different types.
   The ipsi argument works for the S-estimator as well as for the
   MM-estimator.

   - Added implementation of M-S algorithm.

   - Modified subsampling behaviour: avoiding singular resamples by using
     customized LU decomposition.

   - Replaced C style matrices with Fortran style matrices, with as little
     copying as possible.

   - Using LAPACK's DGELS instead of local lu() decomposition.

   - Code clean up: removed all subroutines that were unused.
*/

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <Rmath.h>
#include <complex.h>
// #include <R.h>

#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#ifndef FCONE
# define FCONE
#endif

#include "RobStatTM.h"
//-> <R.h>, <Rinternals.h>  -> XLENGTH, R_xlen_t


/* these  will also move to "lmrob.h" ---
 *  but first make many of these 'static' <<< FIXME!
 */
void fast_s_large_n(double *X, double *y,
		    int *nn, int *pp, int *nRes, int *max_it_scale,
		    int *ggroups, int *nn_group,
		    int *K, int *max_k, double rel_tol, double inv_tol, int *converged,
		    int *best_r, double *bb, double *rrhoc, int *iipsi,
		    double *bbeta, double *sscale, int trace_lev, int mts, int ss);

void fast_s(double *X, double *y,
	    int *nn, int *pp, int *nRes, int *max_it_scale,
	    int *K, int *max_k, double rel_tol, double inv_tol, int *converged,
	    int *best_r, double *bb, double *rrhoc, int *iipsi,
	    double *bbeta, double *sscale, int trace_lev, int mts, int ss);

Rboolean rwls(const double X[], const double y[], int n, int p,
	 double *estimate, double *i_estimate,
	 double *resid, double *loss,
	 double scale, double epsilon,
	 int *max_it, double *rho_c, const int ipsi, int trace_lev);

static void sample_noreplace(int *x, int n, int k, int *ind_space);


double norm2 (double *x, int n);
double norm (double *x, int n);
double norm1(double *x, int n);
double norm_diff2 (double *x, double *y, int n);
double norm_diff (double *x, double *y, int n);
double norm1_diff(double *x, double *y, int n);
double normcnst(const double c[], int ipsi);
double rho_inf (const double c[], int ipsi);

double rho(double x, const double c[], int ipsi);
double psi(double x, const double c[], int ipsi);
double psip(double x, const double c[], int ipsi);// psi'
double psi2(double x, const double c[], int ipsi);// psi''
double wgt(double x, const double c[], int ipsi);

double rho_huber(double x, const double c[]);
double psi_huber(double x, const double c[]);
double psip_huber(double x, const double c[]);
double psi2_huber(double x, const double c[]);
double wgt_huber(double x, const double c[]);

double rho_biwgt(double x, const double c[]);
double psi_biwgt(double x, const double c[]);
double psip_biwgt(double x, const double c[]);
double psi2_biwgt(double x, const double c[]);
double wgt_biwgt(double x, const double c[]);

double rho_gwgt(double x, const double c[]);
double psi_gwgt(double x, const double c[]);
double psip_gwgt(double x, const double c[]);
double wgt_gwgt(double x, const double c[]);

/* PHIONE := phi(1.0) = dnorm(1.0) */
#define PHIONE 0.2419707245191433653275
/* imanginary error function implemented in erfz.c */
double complex erfi(double complex z);
/* antiderivative of optimal psi in support interval */
double Psi_opt(double x, const double c[]);

double rho_opt(double x, const double c[]);
double psi_opt(double x, const double c[]);
double psip_opt(double x, const double c[]);
double wgt_opt(double x, const double c[]);

double rho_optv0(double x, const double c[]);

double rho_hmpl(double x, const double c[]);
double psi_hmpl(double x, const double c[]);
double psip_hmpl(double x, const double c[]);
double psi2_hmpl(double x, const double c[]);
double wgt_hmpl(double x, const double c[]);

double rho_ggw(double x, const double c[]);
void psi_ggw_vec(double *x, int n, void *k);
double psi_ggw(double x, const double c[]);
double psip_ggw(double x, const double c[]);
double wgt_ggw(double x, const double c[]);

double rho_lqq(double x, const double c[]);
double psi_lqq(double x, const double c[]);
double psip_lqq(double x, const double c[]);
double psi2_lqq(double x, const double c[]);
double wgt_lqq(double x, const double c[]);

double rho_modOpt(double x, const double c[]);
double psi_modOpt(double x, const double c[]);
double psip_modOpt(double x, const double c[]);
double wgt_modOpt(double x, const double c[]);

double rho_modOptv0(double x, const double c[]);

double sum_rho_sc(const double r[], double scale, int n, int p,
		  const double c[], int ipsi);
void get_weights_rhop(const double r[], double s, int n, const double rrhoc[], int ipsi,
		      /* --> */ double *w);
int refine_fast_s(const double X[], double *wx, const double y[], double *wy,
		  double *weights, int n, int p, double *res,
		  double *work, int lwork,
		  double *beta_cand,
		  int kk, Rboolean *conv, int max_k, double rel_tol,
		  int trace_lev,
		  double b, double *rrhoc, int ipsi, double initial_scale,
		  double *beta_ref, double *scale);

void m_s_subsample(double *X1, double *y, int n, int p1, int p2,
		   int nResample, int max_it_scale, double rel_tol, double inv_tol, double *bb,
		   double *rrhoc, int ipsi, double *sscale, int trace_lev,
		   double *b1, double *b2, double *t1, double *t2,
		   double *y_tilde, double *res, double *x1, double *x2,
		   int *NIT, int *K, int *KODE, double *SIGMA, double *BET0,
		   double *SC1, double *SC2, double *SC3, double *SC4, int mts, Rboolean ss);

Rboolean m_s_descent(double *X1, double *X2, double *y,
		 int n, int p1, int p2, int K_m_s, int max_k, int max_it_scale,
		 double rel_tol, double *bb, double *rrhoc, int ipsi,
		 double *sscale, int trace_lev,
		 double *b1, double *b2, double *t1, double *t2,
		 double *y_tilde, double *res, double *res2, double *x1, double *x2,
		 int *NIT, int *K, int *KODE, double *SIGMA,  double *BET0,
		 double *SC1, double *SC2, double *SC3, double *SC4);

Rboolean subsample(const double x[], const double y[], int n, int m,
		   double *beta, int *ind_space, int *idc, int *idr,
		   double *lu, double *v, int *p,
		   double *Dr, double *Dc, int rowequ, int colequ,
		   Rboolean sample, int mts, Rboolean ss, double tol_inv, Rboolean solve);

int fast_s_with_memory(double *X, double *y,
		       int *nn, int *pp, int *nRes, int *max_it_scale,
		       int *K, int *max_k, double rel_tol, double inv_tol,  int trace_lev,
		       int *best_r, double *bb, double *rrhoc, int *iipsi,
		       double **best_betas, double *best_scales, int mts, int ss);

/* for "tracing" only : */
void disp_mat(double **a, int n, int m);
void disp_vec(double *a, int n);
void disp_veci(int *a, int n);

double kthplace(double *, int, int);

int find_max(double *a, int n);

double find_scale(double *r, double b, double *rrhoc, int ipsi,
		  double initial_scale, int n, int p, int max_iter);

double median_abs(double *, int, double *);
double MAD(double *a, int n, double center, double *tmp, double *tmp2);

void zero_mat(double **a, int n, int m);


#define INIT_WLS(_X_, _y_, _n_, _p_)                            \
    /* Determine optimal block size for work array*/            \
    F77_CALL(dgels)("N", &_n_, &_p_, &one, _X_, &_n_, _y_,      \
		    &_n_, &work0, &lwork, &info FCONE);               \
    if (info) {                                                 \
	warning(" Problem determining optimal block size, using minimum"); \
	lwork = 2*_p_;                                          \
    } else                                                      \
	lwork = (int)work0;                                     \
                                                                \
    if (trace_lev >= 4)                                         \
	Rprintf(" Optimal block size for DGELS: %d\n", lwork); \
                                                                \
    /* allocate */                                              \
    work =    (double *) Calloc(lwork, double);                 \
    weights = (double *) Calloc(n,     double);

#define CLEANUP_WLS                                             \
    Free(work); Free(weights);

#define CLEANUP_EQUILIBRATION                                   \
    Free(Dr); Free(Dc); Free(Xe);

#define CLEANUP_SUBSAMPLE                                       \
    Free(ind_space); Free(idc); Free(idr); Free(pivot);         \
    Free(lu); Free(v);                                          \
    CLEANUP_EQUILIBRATION;

#define FIT_WLS(_X_, _x_, _y_, _n_, _p_)			\
    /* add weights to _y_ and _x_ */                            \
    for (j=0; j<_n_; j++) {                                     \
    wtmp = sqrt(weights[j]);                                    \
    _y_[j] *= wtmp;                                             \
    for (k=0; k<_p_; k++)                                       \
	_x_[_n_*k+j] = _X_[_n_*k+j] * wtmp;                     \
    }                                                           \
    /* solve weighted least squares problem */                  \
    F77_CALL(dgels)("N", &_n_, &_p_, &one, _x_, &_n_, _y_,      \
		    &_n_, work, &lwork, &info FCONE);                 \
    if (info) {					                \
	if (info < 0) {                                         \
	    CLEANUP_WLS;					\
	    error("DGELS: illegal argument in %i. argument.", info); \
	} else {                                                \
	    if (trace_lev >= 4) {				\
		Rprintf(" Robustness weights in failing step: ");	\
		disp_vec(weights, _n_);				\
	    }                                                   \
	    CLEANUP_WLS;					\
	    error("DGELS: weighted design matrix not of full rank (column %d).\nUse control parameter 'trace.lev = 4' to get diagnostic output.", info); \
	}                                                       \
    }

#define SETUP_EQUILIBRATION(_n_, _p_, _X_, _large_n_)           \
    /* equilibration of matrix _X_                          */  \
    /* solve (Dr X Dc) b = Dr y with beta = Dc b instead of */  \
    /*            X beta = y                                */  \
    /* see Demmel (1997) APPLIED NUMERICAL LINEAR ALGEBRA   */  \
    /*     Section 2.5.2 Equilibration                      */  \
    double *Dr, *Dc, *Xe, rowcnd, colcnd, amax;			\
    int rowequ = 0 , colequ = 0;                                \
    Dr =        (double *) Calloc(_n_,     double);             \
    Dc =        (double *) Calloc(_p_,     double);             \
    Xe =        (double *) Calloc(_n_*_p_, double);             \
    COPY(_X_, Xe, _n_*_p_);                                     \
    F77_CALL(dgeequ)(&_n_, &_p_, Xe, &_n_, Dr, Dc, &rowcnd,	\
    		     &colcnd, &amax, &info);                    \
    if (info) {                                                 \
	if (info < 0) {                                         \
	    CLEANUP_EQUILIBRATION;				\
	    error("DGEEQ: illegal argument in %i. argument", -1 * info); \
	} else if (info > _n_) {                                \
	    if (_large_n_) {                                    \
	        error("Fast S large n strategy failed. Use control parameter 'fast.s.large.n = Inf'."); \
	    } else {						\
                error("DGEEQU: column %i of the design matrix is exactly zero.", info - _n_); \
	    }                                                   \
	} else {                                                \
	/* FIXME: replace dgeequ by our own version */          \
	/* that does not treat this as error */                 \
	    warning(" Skipping design matrix equilibration (DGEEQU): row %i is exactly zero.", info); \
	}                                                       \
    } else {							\
        /* scale _X_ */                                         \
        char equed = ' ';         					\
	F77_CALL(dlaqge)(&_n_, &_p_, Xe, &_n_, Dr, Dc, &rowcnd,	\
			 &colcnd, &amax, &equed FCONE);		\
        rowequ = equed == 'B' || equed == 'R';                  \
	colequ = equed == 'B' || equed == 'C';                  \
    }

#define SETUP_SUBSAMPLE(_n_, _p_, _X_, _large_n_)		\
    /* (Pointers to) Arrays - to be allocated */                \
    int *ind_space, *idc, *idr, *pivot;				\
    double *lu, *v;						\
    ind_space = (int *)    Calloc(_n_,     int);                \
    idc =       (int *)    Calloc(_n_,     int);                \
    idr =       (int *)    Calloc(_p_,     int);                \
    pivot =     (int *)    Calloc(_p_-1,   int);                \
    lu =        (double *) Calloc(_p_*_p_, double);             \
    v =         (double *) Calloc(_p_,     double);             \
    SETUP_EQUILIBRATION(_n_, _p_, _X_, _large_n_);

#define COPY(from, to, len) Memcpy(to, from, len)
/* This assumes that 'p' is correctly defined, and 'j' can be used in caller: */
/* #define COPY(BETA_FROM, BETA_TO, _p_)			\ */
/*     for(j=0; j < _p_; j++) BETA_TO[j] = BETA_FROM[j]; */
/* In theory BLAS should be fast, but this seems slightly slower,
 * particularly for non-optimized BLAS :*/
/* static int one = 1; */
/* #define COPY(BETA_FROM, BETA_TO, _p_) \ */
/*     F77_CALL(dcopy)(&_p_, BETA_FROM, &one, BETA_TO, &one);  */

#define EPS_SCALE 1e-10
#define INFI 1e+20

/* Called from R, this function computes an S-regression estimator */
void R_lmrob_S(double *X, double *y, int *n, int *P,
	       int *nRes, // = nResample ( = 500, by default)
	       double *scale, double *beta_s,
	       double *rrhoc, int *iipsi, double *bb,
	       int *best_r, int *Groups, int *N_group,
	       int *K_s, int *max_k, int *max_it_scale, //double *rel_tol_scale,
	       double *rel_tol, double *inv_tol, int *converged,
	       int *trace_lev, int *mts, int *ss, int *cutoff)
{
    /* best_r = 't' of Salibian-Barrera_Yohai(2006),
     *	      = no. of best candidates to be iterated further ("refined")
     *        =  2, by default
     */

    if (*nRes > 0) {
	if (*n > *cutoff) {
	    if(*trace_lev > 0)
		Rprintf("lmrob_S(n = %d, nRes = %d): fast_s_large_n():\n", *n, *nRes);
	    fast_s_large_n(X, y, n, P, nRes, max_it_scale,
			   Groups, N_group,
			   K_s, max_k, *rel_tol, *inv_tol, converged,
			   best_r, bb, rrhoc, iipsi, beta_s, scale, *trace_lev, *mts, *ss);
	} else {
	    if(*trace_lev > 0)
		Rprintf("lmrob_S(n = %d, nRes = %d): fast_s() [non-large n]:\n", *n, *nRes);
	    fast_s(X, y, n, P, nRes, max_it_scale,
		   K_s, max_k, *rel_tol, *inv_tol, converged,
		   best_r, bb, rrhoc, iipsi, beta_s, scale, *trace_lev, *mts, *ss);
	}
    } else {
	if(*trace_lev > 0)
	    Rprintf("lmrob_S(nRes = 0, n = %d): --> find_scale() only:\n", *n);
	*scale = find_scale(y, *bb, rrhoc, *iipsi, *scale, *n, *P,
			    *max_it_scale);
    }
}

/* Called from R, this function computes an M-S-regression estimator */
// not only called from ../R/lmrob.M.S.R,  but also ../inst/xtraR/m-s_fns.R
//			~~~~~~~~~~~~~~~~	    ~~~~~~~~~~~~~~~~~~~~~~~
void R_lmrob_M_S(double *X1, double *X2, double *y, double *res,
		 int *nn, int *pp1, int *pp2, int *nRes, int *max_it_scale,
		 double *scale, double *b1, double *b2,
		 double *rho_c, int *ipsi, double *bb,
		 int *K_m_s, int *max_k, double *rel_tol, double *inv_tol,
		 int *converged,
		 int *trace_lev,
		 int *orthogonalize, int *subsample, int *descent,
		 int *mts, int *ss)
{
    /* Initialize (some of the) memory here,
     * so that we have to do it only once */
    int i, n = *nn, p1 = *pp1, p2 = *pp2, one = 1;

    /* (Pointers to) Arrays - to be allocated */
    double *t1, *t2, *y_tilde, *y_work, done = 1., dmone = -1.;
    double *x1, *x2, *ot1, *oT2, *ptr;

    if(*trace_lev > 0) Rprintf(
	"lmrob_M_S(n = %d, nRes = %d, (p1,p2)=(%d,%d), (orth,subs,desc)=(%d,%d,%d))\n",
	n, *nRes, p1, p2,
	*orthogonalize, *subsample, *descent);

    t1 =      (double *) R_alloc(n,  sizeof(double)); /* size n needed for rllarsbi */
    t2 =      (double *) R_alloc(p2, sizeof(double));
    ot1 =     (double *) R_alloc(p1, sizeof(double));
    oT2 =     (double *) R_alloc(p2*p1, sizeof(double));
    y_work =  (double *) R_alloc(n,  sizeof(double));
    COPY(y, y_work, n);
    y_tilde = (double *) R_alloc(n,  sizeof(double));
    x1 =      (double *) R_alloc(n*p1, sizeof(double));
    x2 =      (double *) R_alloc(n*p2, sizeof(double));
    COPY(X2, x2, n*p2);

    /* Variables required for rllarsbi
     *  (l1 / least absolut residuals - estimate) */
    int NIT=0, K=0, KODE=0;
    double SIGMA = 0.,
	*SC1 = (double *) R_alloc(n,  sizeof(double)),
	*SC2 = (double *) R_alloc(p1, sizeof(double)),
	*SC3 = (double *) R_alloc(p1, sizeof(double)),
	*SC4 = (double *) R_alloc(p1, sizeof(double));
    double BET0 = 0.773372647623; /* = pnorm(0.75) */

    /* STEP 1: Orthgonalize X2 and y from X1 */
    if (*orthogonalize) {
	COPY(X1, x1, n*p1);
	F77_CALL(rllarsbi)(x1, y_work, &n, &p1, &n, &n, rel_tol,
			   &NIT, &K, &KODE, &SIGMA, t1, y_tilde, SC1, SC2,
			   SC3, SC4, &BET0);
	COPY(t1, ot1, p1);
	for (i=0; i < p2; i++) {
	    COPY(X1, x1, n*p1);
	    ptr = X2+i*n; COPY(ptr, y_work, n);
	    F77_CALL(rllarsbi)(x1, y_work, &n, &p1, &n, &n, rel_tol,
			       &NIT, &K, &KODE, &SIGMA, t1, x2+i*n, SC1, SC2,
			       SC3, SC4, &BET0);
	    ptr = oT2+i*p1; COPY(t1, ptr, p1);
	}
	COPY(y_tilde, y_work, n);
	/* compare with Maronna & Yohai 2000:
	 * y_work and y_tilde now contain \tilde y, ot1 -> t_1,
	 * x2 -> \tilde x2, oT2 -> T_2 */
    }

    /* STEP 2: Subsample */
    if (*subsample) {
	m_s_subsample(X1, y_work, n, p1, p2, *nRes, *max_it_scale,
		      *rel_tol, *inv_tol, bb,
		      rho_c, *ipsi, scale, *trace_lev,
		      b1, b2, t1, t2, y_tilde, res, x1, x2,
		      &NIT, &K, &KODE, &SIGMA, &BET0,
		      SC1, SC2, SC3, SC4, *mts, *ss);

	if (*scale < 0)
	    error("m_s_subsample() stopped prematurely (scale < 0).");
    }

    /* STEP 3: Transform back */
    if (*orthogonalize) {
	/* t1 = ot1 + b1 - oT2 %*% b2 */
	for(int i=0; i < p1; i++) t1[i] = ot1[i] + b1[i];
	F77_CALL(dgemv)("N", &p1, &p2, &dmone, oT2, &p1, b2, &one, &done, t1, &one FCONE);
	COPY(t1, b1, p1);
	/* restore x2 */
	COPY(X2, x2, n*p2);
    }

    /* update / calculate residuals */
    COPY(y, res, n);
    F77_CALL(dgemv)("N", &n, &p1, &dmone, X1, &n, b1, &one, &done, res, &one FCONE);
    F77_CALL(dgemv)("N", &n, &p2, &dmone, X2, &n, b2, &one, &done, res, &one FCONE);

    /* STEP 4: Descent procedure */
    if (*descent) {
	*converged = m_s_descent(
	    X1, X2, y, n, p1, p2, *K_m_s, *max_k, *max_it_scale,
	    *rel_tol, bb, rho_c, *ipsi, scale, *trace_lev,
	    b1, b2, t1, t2, y_tilde, res, y_work, x1, x2,
	    &NIT, &K, &KODE, &SIGMA, &BET0, SC1, SC2, SC3, SC4);
    }
}

/* This function performs RWLS iterations starting from
 * an S-regression estimator (and associated residual scale).
 * So, in itself, this is ``just'' an M-estimator -- called from R's
 * lmrob..M..fit()  [ ../R/lmrob.MM.R ]
 * ~~~~~~~~~~~~~~~
 * NOTE: rel_tol now controls the *relative* changes in beta,
 *      instead of being hard-wired to EPS = 1e-7 and bounding the
 *	absolute || beta_1 - beta_2 ||
 */
void R_lmrob_MM(double *X, double *y, int *n, int *P,
		double *beta_initial, double *scale,
		double *beta_m, double *resid,
		int *max_it, double *rho_c, int *ipsi, double *loss,
		double *rel_tol, int *converged, int *trace_lev, int *mts, int *ss)
{
    /* starting from the S-estimate (beta_initial), use
     * irwls to compute the MM-estimate (beta_m)  */

    if(*trace_lev > 0)
	Rprintf("lmrob_MM(): rwls():\n");

    *converged = (int)rwls(X,y,*n,*P,beta_m, beta_initial, resid, loss,
		      *scale, *rel_tol,
		      max_it, rho_c, *ipsi, *trace_lev);
    if (!converged)
	COPY(beta_initial, beta_m, *P);
}


/* Call subsample() from R, for testing purposes only */
void R_subsample(const double x[], const double y[], int *n, int *m,
		 double *beta, int *ind_space, int *idc, int *idr,
		 double *lu, double *v, int *p,
		 double *_Dr, double *_Dc, int *_rowequ, int *_colequ,
		 int *status, int *sample, int *mts, int *ss, double *tol_inv,
		 int *solve)
{
    int info;

    /*	set the seed */
    GetRNGstate();

    SETUP_EQUILIBRATION(*n, *m, x, 0);

    *status = subsample(Xe, y, *n, *m, beta, ind_space, idc, idr, lu, v, p,
			Dr, Dc, rowequ, colequ, (Rboolean)*sample, *mts, (Rboolean)*ss,
			*tol_inv, (Rboolean)*solve);

    COPY(Dr, _Dr, *n);
    COPY(Dc, _Dc, *m);
    *_rowequ = rowequ;
    *_colequ = colequ;

    CLEANUP_EQUILIBRATION;

    PutRNGstate();
}

/*----------------------------------------------------------------------------*/

SEXP R_psifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_) {
    /*
     * Calculate psi for vectorized x, scaled to get  psi'(0) = 1
     * deriv -1: rho(x)  {*not* normalized}
     * deriv  0: psi(x)  = rho'(x)
     * deriv  1: psi'(x) = rho''(x)   {we always have psip(0) == 1}
     * deriv  2: psi''(x)= rho'''(x)
     */
    int nprot = 1, ipsi = asInteger(ipsi_), deriv = asInteger(deriv_);
    if (isInteger(x_)) {
	x_ = PROTECT(coerceVector(x_, REALSXP)); nprot++;
    }
    if (!isReal(x_)) error(_("Argument '%s' must be numeric or integer"), "x");
    if (!isReal(c_)) error(_("Argument '%s' must be numeric or integer"), "cc");
    R_xlen_t i, n = XLENGTH(x_);
    SEXP res = PROTECT(allocVector(REALSXP, n)); // the result
    double *x = REAL(x_), *r = REAL(res), *cc = REAL(c_);

// put the for() loop *inside* the switch (--> speed for llength >> 1) :
#define for_i_n_NA  for(i = 0; i < n; i++) r[i] = ISNAN(x[i]) ? x[i] :

    switch(deriv) { // our rho() is rho~(), i.e., scaled to max = 1
    case -1: {
	double rho_Inf = rho_inf(cc, ipsi);
	for_i_n_NA rho(x[i], cc, ipsi) * rho_Inf; break;
    }
    case  0: for_i_n_NA psi (x[i], cc, ipsi); break;
    case  1: for_i_n_NA psip(x[i], cc, ipsi); break;
    case  2: for_i_n_NA psi2(x[i], cc, ipsi); break;
    default: error(_("'deriv'=%d is invalid"), deriv);
    }
    UNPROTECT(nprot);
    return res;
}

SEXP R_chifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_) {
    /*
     * Calculate chi for vectorized x, i.e. rho~(.) with rho~(inf) = 1:
     * deriv 0: chi (x)  = \rho(x) / \rho(Inf) =: \rho(x) * nc  == our rho() C-function
     * deriv 1: chi'(x)  = psi(x)  * nc
     * deriv 2: chi''(x) = psi'(x) * nc
     */
    int nprot = 1, ipsi = asInteger(ipsi_), deriv = asInteger(deriv_);
    if (isInteger(x_)) {
	x_ = PROTECT(coerceVector(x_, REALSXP)); nprot++;
    }
    if (!isReal(x_)) error(_("Argument '%s' must be numeric or integer"), "x");
    if (!isReal(c_)) error(_("Argument '%s' must be numeric or integer"), "cc");
    R_xlen_t i, n = XLENGTH(x_);
    SEXP res = PROTECT(allocVector(REALSXP, n)); // the result
    double *x = REAL(x_), *r = REAL(res), *cc = REAL(c_);

    // our rho() is rho~() == chi(), i.e., scaled to max = 1
    double rI = (deriv > 0) ? rho_inf(cc, ipsi) : 0./* -Wall */;

    switch(deriv) {
    case 0: for_i_n_NA rho(x[i], cc, ipsi); break;
    case 1: for_i_n_NA psi (x[i], cc, ipsi) / rI; break;
    case 2: for_i_n_NA psip(x[i], cc, ipsi) / rI; break;
    case 3: for_i_n_NA psi2(x[i], cc, ipsi) / rI; break;
    default: error(_("'deriv'=%d is invalid"), deriv);
    }
    UNPROTECT(nprot);
    return res;
}

SEXP R_wgtfun(SEXP x_, SEXP c_, SEXP ipsi_) {
    /*
     * Calculate wgt(x) = psi(x)/x for vectorized x
     */
    int nprot = 1, ipsi = asInteger(ipsi_);
    if (isInteger(x_)) {
	x_ = PROTECT(coerceVector(x_, REALSXP)); nprot++;
    }
    if (!isReal(x_)) error(_("Argument '%s' must be numeric or integer"), "x");
    if (!isReal(c_)) error(_("Argument '%s' must be numeric or integer"), "cc");
    R_xlen_t i, n = XLENGTH(x_);
    SEXP res = PROTECT(allocVector(REALSXP, n)); // the result
    double *x = REAL(x_), *r = REAL(res), *cc = REAL(c_);

    for_i_n_NA wgt(x[i], cc, ipsi);

    UNPROTECT(nprot);
    return res;
}

#undef for_i_n_NA

SEXP R_rho_inf(SEXP cc, SEXP ipsi) {
    if (!isReal(cc)) error(_("Argument 'cc' must be numeric"));
    if (!isInteger(ipsi)) error(_("Argument 'ipsi' must be integer"));
    return ScalarReal(rho_inf(REAL(cc), INTEGER(ipsi)[0]));
}

double rho_inf(const double k[], int ipsi) {
    /*
     * Compute  \rho(\infty) for psi functions
     * (Note that our C function rho() is "rho~" and has rho(Inf) = 1)
     */
    double c = k[0];

    switch(ipsi) {
    default: error("rho_inf(): ipsi=%d not implemented.", ipsi);
    case 0: return(R_PosInf); // huber
    case 1: return(c*c/6.); // biweight
    case 2: return(c*c); // GaussWeight / "Welsh"
    case 3: return(PHIONE / (PHIONE - k[0]) * k[3] * k[3] * k[5]); // Optimal poly (same as optimal original)
    case 4: return(0.5*k[0]*(k[1]+k[2]-k[0])); // Hampel
    case 5: // GGW (Generalized Gauss Weight)
	switch((int)c) {
	default:
	case 0: return(k[4]); break; // k[4] == cc[5] in R -- must be correct!
	case 1: return(5.309853); break;
	case 2: return(2.804693); break;
	case 3: return(0.3748076); break;
	case 4: return(4.779906); break;
	case 5: return(2.446574); break;
	case 6: return(0.4007054); break;
	};
    case 6: // LQQ aka 'lin psip'
	return (k[2]*k[1]*(3*k[1]+2*k[0]) + (k[0]+k[1])*(k[0]+k[1])) / (6.*(k[2]-1.)); break;
    case 7: return(k[3] * k[3] * k[5]); break; // Modified optimal poly
    case 8: return(PHIONE / (PHIONE - k[0]) * k[3] * k[3] * k[5]); // Optimalv0
    case 9: return(k[3] * k[3] * k[5]); break; // Modified Optimalv0
      
    }
} // rho_inf()

double normcnst(const double k[], int ipsi) {
    /*
     * return normalizing constant for psi functions
     */

    double c = k[0];

    switch(ipsi) {
    default: error("normcnst(): ipsi=%d not implemented.", ipsi);
    case 0: return(0.); // huber {normcnst() should never be used for that!}
    case 1: return(6./(c*c)); // biweight
    case 2: return(1./(c*c)); // GaussWeight / "Welsh"
    case 3: return(1./3.25/(c*c)); // Optimal FIXME: maybe return NA?
    case 4: return(2./(k[0]*(k[1]+k[2]-k[0]))); // Hampel
    case 5: // GGW
	switch((int)c) {
	default:
	case 0: return(1./ k[4]); break; // k[4] == cc[5] in R -- must be correct!
	case 1: return(1./5.309853); break;
	case 2: return(1./2.804693); break;
	case 3: return(1./0.3748076); break;
	case 4: return(1./4.779906); break;
	case 5: return(1./2.446574); break;
	case 6: return(1./0.4007054); break;
	};
    case 6: // LQQ aka 'lin psip'
	return((6*(k[2]-1))/(k[2]*k[1]*(3*k[1]+2*k[0])+(k[0]+k[1])*(k[0]+k[1])));
    }
} // normcnst()


double rho(double x, const double c[], int ipsi)
{
    /*
     * return the correct rho according to ipsi
     * This rho() is normalized to 1, called rho~() or chi() in other contexts
     */
    switch(ipsi) {
    default: error("rho(): ipsi=%d not implemented.", ipsi);
    case 0: return(rho_huber(x, c)); // huber
    case 1: return(rho_biwgt(x, c)); // biweight
    case 2: return(rho_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(rho_opt(x, c)); // Optimal poly
    case 4: return(rho_hmpl(x, c)); // Hampel
    case 5: return(rho_ggw(x, c)); // GGW (Generalized Gauss Weight)
    case 6: return(rho_lqq(x, c)); // LQQ := Linear-Quadratic-Quadratic
	// was LGW := "lin psip" := piecewise linear psi'()
    case 7: return(rho_modOpt(x, c)); // Modified Optimal poly
    case 8: return(rho_optv0(x, c)); // Optimalv0
    case 9: return(rho_modOptv0(x, c)); // Modified Optimalv0
    }
}

double psi(double x, const double c[], int ipsi)
{
    /*
     * return the correct psi according to ipsi
     * this is actually rho' and not psi
     */
    switch(ipsi) {
    default: error("psi(): ipsi=%d not implemented.", ipsi);
    case 0: return(psi_huber(x, c)); // huber
    case 1: return(psi_biwgt(x, c)); // biweight
    case 2: return(psi_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(psi_opt(x, c)); // Optimal poly same as optimal
    case 4: return(psi_hmpl(x, c)); // Hampel
    case 5: return(psi_ggw(x, c)); // GGW
    case 6: return(psi_lqq(x, c)); // LQQ (piecewise linear psi')
    case 7: return(psi_modOpt(x, c)); // Modified Optimal poly same as modopt
    case 8: return(psi_opt(x, c)); // Optimalv0
    case 9: return(psi_modOpt(x, c)); // Modified Optimalv0
    }
}

double psip(double x, const double c[], int ipsi)
{
    /*
     * return the correct ppsi according to ipsi
     * this is actually rho'' and not psip
     */
    switch(ipsi) {
    default: error("psip(): ipsi=%d not implemented.", ipsi);
    case 0: return(psip_huber(x, c)); // huber
    case 1: return(psip_biwgt(x, c)); // biweight
    case 2: return(psip_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(psip_opt(x, c)); // Optimal poly same as optimal
    case 4: return(psip_hmpl(x, c)); // Hampel
    case 5: return(psip_ggw(x, c)); // GGW
    case 6: return(psip_lqq(x, c)); // LQQ (piecewise linear psi')
    case 7: return(psip_modOpt(x, c)); // Modified Optimal poly same as modopt
    case 8: return(psip_opt(x, c)); // Optimalv0
    case 9: return(psip_modOpt(x, c)); // Modified Optimalv0
    }
}

double psi2(double x, const double c[], int ipsi)
{
    /* Compute   psi''(x) == rho'''(x)
     */
    switch(ipsi) {
    // default: error("psi2: ipsi=%d not implemented.", ipsi);
    case 0: return(psi2_huber(x, c)); // huber
    case 1: return(psi2_biwgt(x, c)); // biweight
    case 4: return(psi2_hmpl(x, c)); // Hampel
    case 6: return(psi2_lqq(x, c)); // LQQ (piecewise linear psi')

    default: error("psi2(): ipsi=%d not implemented.", ipsi);
/*
    case 2: return(psi2_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(psi2_opt(x, c)); // Optimal
    case 5: return(psi2_ggw(x, c)); // GGW
    case 7: return(psi2_modOpt(x, c)); // Modified Optimal
*/
    }
}

double wgt(double x, const double c[], int ipsi)
{
    /*
     * return the correct wgt according to ipsi
     * wgt: rho'(x) / x
     */
    switch(ipsi) {
    default:
    case 0: return(wgt_huber(x, c)); // huber
    case 1: return(wgt_biwgt(x, c)); // biweight
    case 2: return(wgt_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(wgt_opt(x, c)); // Optimal  poly same as optimal
    case 4: return(wgt_hmpl(x, c)); // Hampel
    case 5: return(wgt_ggw(x, c)); // GGW
    case 6: return(wgt_lqq(x, c)); // LQQ (piecewise linear psi')
    case 7: return(wgt_modOpt(x, c)); // Modified Optimal poly same as modopt
    case 8: return(wgt_opt(x, c)); // Optimal 
    case 9: return(wgt_modOpt(x, c)); // Modified Optimal
      
    }
}

//---  Huber's rho / psi / ...
//---  -------

/* Huber's rho():  contrary to all the redescenders below,
   this can NOT be scaled to rho(Inf)=1 : */
double rho_huber(double x, const double c[])
{
    return (fabs(x) <= c[0]) ? x*x*0.5 : c[0]*(fabs(x) - c[0]/2);
}

double psi_huber(double x, const double c[])
{
// Huber's psi = rho'()
    return (x <= -c[0]) ? -c[0] : ((x < c[0]) ? x : c[0]);
}

double psip_huber(double x, const double c[])
{
// psi' = rho'' : Second derivative of Huber's loss function
    return (fabs(x) >= c[0]) ? 0. : 1.;
}

double psi2_huber(double x, const double c[])
{
// psi'' = rho''' : Third derivative of Huber's loss function
    return 0;
// FIXME? return NaN when  |x| == c ?? -- then also for psi2_hmpl()
}

double wgt_huber(double x, const double c[])
{
/*
 * Weights for Huber's loss function w(x) = psi(x)/x
 */
    return (fabs(x) >= c[0]) ? c[0]/fabs(x) : 1.;
}

//--- Biweight = Bisquare = Tukey's Biweight ...
//--- --------------------------------------

double rho_biwgt(double x, const double c[])
{
/*
 * Tukey's bisquare loss function  == R's  tukeyChi()
 */
    if (fabs(x) > (*c))
	return(1.);
    else {
	double t = x / (*c);
	t *= t; /* = t^2 */
	return( t*(3. + t*(-3. + t)) );
    }
}

double psi_biwgt(double x, const double c[])
{
/*
 * First derivative of Tukey's bisquare loss function
 */
    if (fabs(x) > (*c))
	return(0.);
    else {
	double a = x / (*c),
	    u = 1. - a*a;
	return( x * u * u );
    }
}

double psip_biwgt(double x, const double c[])
{
/*
 * Second derivative of Tukey's bisquare loss function
 */
    if (fabs(x) > (*c))
	return(0.);
    else {
	x /= *c;
	double x2 = x*x;
	return( (1. - x2) * (1 - 5 * x2));
    }
}

double psi2_biwgt(double x, const double c[])
{
/** 3rd derivative of Tukey's bisquare loss function rho()
 *= 2nd derivative of psi() :
 */
    if (fabs(x) >= c[0]) // psi''()  *is* discontinuous at x = c[0]: use "middle" value there:
	return (fabs(x) == c[0]) ? 4*x/c[0] : 0.;
    else {
	x /= c[0];
	double x2 = x*x;
	return 4*x/c[0] * (5 * x2 - 3.);
    }
}

double wgt_biwgt(double x, const double c[])
{
/*
 * Weights for Tukey's bisquare loss function
 */
    if( fabs(x) > *c )
	return(0.);
    else {
	double a = x / (*c);
	a = (1. - a)*(1. + a);
	return( a * a );
    }
}

//---------- gwgt == Gauss Weight Loss function =: "Welsh" --------------------

double rho_gwgt(double x, const double c[])
{
    /*
     * Gauss Weight Loss function
     */
    double ac = x / (*c);
    return(-expm1(-(ac*ac)/2));
}

// Largest x  such that  exp(-x) does not underflow :
static double MIN_Exp = -708.4; // ~ = M_LN2 * DBL_MIN_EXP = -log(2) * 1022 = -708.3964 */

// Largest x  such that  exp(-x^2/2) does not underflow :
static double MAX_Ex2 = 37.7; // ~ = sqrt(- 2. * M_LN2 * DBL_MIN_EXP);
/* max {x | exp(-x^2/2) < .Machine$double.xmin } =
 * min {x |  x^2 > -2*log(2)* .Machine$double.min.exp } =
 *  = sqrt(-2*log(2)* .Machine$double.min.exp) = {IEEE double}
 *  = sqrt(log(2) * 2044) = 37.64031 */

double psi_gwgt(double x, const double c[])
{
    /*
     * Gauss Weight Psi()
     */
    double a = x / (*c);
    if(fabs(a) > MAX_Ex2)
	return 0.;
    else
	return x*exp(-(a*a)/2);
}

double psip_gwgt(double x, const double c[])
{
    /*
     * Gauss Weight Psi'()
     */
    x /= (*c);
    if(fabs(x) > MAX_Ex2)
	return 0.;
    else {
	double ac = x*x;
	return exp(-ac/2) * (1. - ac);
    }
}

double wgt_gwgt(double x, const double c[])
{
    /*
     * Gauss Weight Loss function
     */
    double a = x / (*c);
    return(exp(-(a*a)/2));
}

/* antiderivative of optimal psi in support interval */
double Psi_opt(double x, const double c[])
{
  return(0.5 * x * x - c[0] * M_PI * creal(erfi((double complex) (M_SQRT1_2 * x))));
}

double rho_optv0(double x, const double c[])
{
  /*
   * c[] contains 6 values:
   *   c[0]: tuning constant a
   *   c[1] & c[2]: endpoints of support interval calculated from a (x > 0)
   *   c[3]: rescaling parameter; 1.0 gives definition of optimal psi
   *   c[4]: Psi_opt(c[1], c)
   *   c[5]: rho_opt(Inf) when c[3] = 1
   */
  
  x = fabs(x) / c[3];
  
  if(x <= c[1])
    return(0.0);
  
  if(x >= c[2])
    return(1.0);
  
  return((Psi_opt(x, c) - c[4]) / c[5]);
}




double rho_opt(double x, const double c[])
{
  /*
   * c[] contains 16 values:
   *   c[0]: tuning constant a
   *   c[1] & c[2]: endpoints of support interval calculated from a (x > 0)
   *   c[3]: rescaling parameter; 1.0 gives definition of optimal psi
   *   c[4]: Psi_opt(c[1], c)
   *   c[5]: rho_opt(Inf) when c[3] = 1
   *   c[6:10]: polynomial coefficients
   *   c[11:12]: endpoints of support (for x > 0)
   *   c[13]: rescaling parameter, 1.0 yields poly optimal psi
   *   c[14], c[15]: reference points (max value is c[15] - c[14])
   */
  
  x = fabs(x) / c[13];
  
  if(x <= c[11])
    return(0.0);
  
  if(x >= c[12])
    return(c[15]-c[14]);
  
  double x2 = x * x;
  
  double u3 = c[6] * x + c[7] * x2 / 2.0 + c[8] * x2 * x2 / 4.0 + c[9] * x2 * x2 * x2 / 6.0 + c[10] * x2 * x2 * x2 * x2 / 8.0;
  
  return(u3 - c[14]);  
}


double psi_opt(double x, const double c[])
{
 /*
  * c[] contains 6 values:
  *   c[0]: tuning constant a
  *   c[1] & c[2]: endpoints of support interval calculated from a (x > 0)
  *   c[3]: rescaling parameter; 1.0 gives definition of optimal psi
  *   c[4]: Psi_opt(c[1], c)
  *   c[5]: rho_opt(Inf) when c[3] = 1
  */

  double absx = fabs(x) / c[3], y = 0.0;

  if(absx > c[1] && absx < c[2]) {
    y = c[3] * PHIONE / (PHIONE - c[0]) * (absx - c[0] / dnorm(absx, 0.0, 1.0, 0));
    y = x > 0.0 ? y : -1.0 * y;
  }

  return(y);
}

double psip_opt(double x, const double c[])
{
 /*
  * c[] contains 6 values:
  *   c[0]: tuning constant a
  *   c[1] & c[2]: endpoints of support interval calculated from a (x > 0)
  *   c[3]: rescaling parameter; 1.0 gives definition of optimal psi
  *   c[4]: Psi_opt(c[1], c)
  *   c[5]: rho_opt(Inf) when c[3] = 1
  */

  x = fabs(x) / c[3];

  if(x > c[1] && x < c[2])
    return(PHIONE / (PHIONE - c[0]) * (1.0 - c[0] * x / dnorm(x, 0.0, 1.0, 0)));

  return(0.0);
}

double wgt_opt(double x, const double c[])
{
 /*
  * c[] contains 6 values:
  *   c[0]: tuning constant a
  *   c[1] & c[2]: endpoints of support interval calculated from a (x > 0)
  *   c[3]: rescaling parameter; 1.0 gives definition of optimal psi
  *   c[4]: Psi_opt(c[1], c)
  *   c[5]: rho_opt(Inf) when c[3] = 1
  */

  x = fabs(x) / c[3];

  if(x > c[1] && x < c[2])
    return(PHIONE / (PHIONE - c[0]) * (1.0 - c[0] / (x * dnorm(x, 0.0, 1.0, 0))));

  return(0.0);
}


double rho_hmpl(double x, const double k[])
{
    /*
     * rho()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     *
     * This function is normalized s.t. rho(inf) = 1
     */
    double u = fabs(x),
	nc = k[0]*(k[1]+k[2]-k[0])/2;

    if (u <= k[0])
	return( x*x/2 / nc );
    else if (u <= k[1])
	return( ( u - k[0]/2 ) * k[0] / nc );
    else if (u <= k[2])
	return( ( k[1] - k[0]/2 + (u - k[1]) * (1 - ( u - k[1] ) / ( k[2] - k[1] ) / 2 )) * k[0] / nc);
    else
	return( 1 );
}

double psi_hmpl(double x, const double k[])
{
    /*
     * psi()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    // double sx = sign(x), u = fabs(x); :
    double sx, u;
    if (x < 0) { sx = -1; u = -x; } else { sx = +1; u = x; }

    if (u <= k[0])
	return( x );
    else if (u <= k[1])
	return sx * k[0];
    else if (u <= k[2])
	return sx * k[0] * (k[2] - u) / (k[2] - k[1]);
    else
	return 0.;
}

double psip_hmpl(double x, const double k[])
{
    /*
     * psi'()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    double u = fabs(x);

    if (u <= k[0])
	return( 1 );
    else if (u <= k[1])
	return( 0 );
    else if (u <= k[2])
	return( k[0] / ( k[1] - k[2] ) );
    else
	return( 0 );
}

double psi2_hmpl(double x, const double k[])
{
    /*
     * psi''()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    return 0.; // even though psi'() is already discontinuous at k[j]
}

double wgt_hmpl(double x, const double k[])
{
    /*
     * w(x) = psi(x)/x  for Hampel's redescending psi function
     * Hampel redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    double u = fabs(x);

    if (u <= k[0])
	return( 1 );
    else if (u <= k[1])
	return( k[0] / u );
    else if (u <= k[2])
	return( k[0] * ( k[2] - u ) / ( k[2] - k[1] ) / u );
    else
	return( 0 );
}


//--- GGW := Generalized Gauss-Weight    Koller and Stahel (2011)
//--- ---
// rho() & chi()  need to be calculated by numerical integration

double rho_ggw(double x, const double k[])
{
    /*
     * Gauss Weight with constant center
     */
    if (k[0] > 0) { // for hard-coded constants
	const double C[6][20] = { // 0: b = 1, 95% efficiency
	    {0.094164571656733, -0.168937372816728, 0.00427612218326869,
	     0.336876420549802, -0.166472338873754, 0.0436904383670537,
	     -0.00732077121233756, 0.000792550423837942, -5.08385693557726e-05,
	     1.46908724988936e-06, -0.837547853001024, 0.876392734183528,
	     -0.184600387321924, 0.0219985685280105, -0.00156403138825785,
	     6.16243137719362e-05, -7.478979895101e-07, -3.99563057938975e-08,
	     1.78125589532002e-09, -2.22317669250326e-11},
	    // 1: b = 1, 85% efficiency
	    {0.174505224068561, -0.168853188892986, 0.00579250806463694,
	     0.624193375180937, -0.419882092234336, 0.150011303015251,
	     -0.0342185249354937, 0.00504325944243195, -0.0004404209084091,
	     1.73268448820236e-05, -0.842160072154898, 1.19912623576069,
	     -0.345595777445623, 0.0566407000764478, -0.00560501531439071,
	     0.000319084704541442, -7.4279004383686e-06, -2.02063746721802e-07,
	     1.65716101809839e-08, -2.97536178313245e-10},
	    // 2: b = 1, bp 0.5
	    {1.41117142330711, -0.168853741371095, 0.0164713906344165,
	     5.04767833986545, -9.65574752971554, 9.80999125035463,
	     -6.36344090274658, 2.667031271863, -0.662324374141645,
	     0.0740982983873332, -0.84794906554363, 3.4315790970352,
	     -2.82958670601597, 1.33442885893807, -0.384812004961396,
	     0.0661359078129487, -0.00557221619221031, -5.42574872792348e-05,
	     4.92564168111658e-05, -2.80432020951381e-06},
	    // 3: b = 1.5, 95% efficiency
	    {0.104604570079252, 0.0626649856211545, -0.220058184826331,
	     0.403388189975896, -0.213020713708997, 0.102623342948069,
	     -0.0392618698058543, 0.00937878752829234, -0.00122303709506374,
	     6.70669880352453e-05, 0.632651530179424, -1.14744323908043,
	     0.981941598165897, -0.341211275272191, 0.0671272892644464,
	     -0.00826237596187364, 0.0006529134641922, -3.23468516804340e-05,
	     9.17904701930209e-07, -1.14119059405971e-08},
	    // 4: b = 1.5, 85% efficiency
	    {0.205026436642222, 0.0627464477520301, -0.308483319391091,
	     0.791480474953874, -0.585521414631968, 0.394979618040607,
	     -0.211512515412973, 0.0707208739858416, -0.0129092527174621,
	     0.000990938134086886, 0.629919019245325, -1.60049136444912,
	     1.91903069049618, -0.933285960363159, 0.256861783311473,
	     -0.0442133943831343, 0.00488402902512139, -0.000338084604725483,
	     1.33974565571893e-05, -2.32450916247553e-07},
	    // 5: b = 1.5, bp 0.5
	    {1.35010856132000, 0.0627465630782482, -0.791613168488525,
	     5.21196700244212, -9.89433796586115, 17.1277266427962,
	     -23.5364159883776, 20.1943966645350, -9.4593988142692,
	     1.86332355622445, 0.62986381140768, -4.10676399816156,
	     12.6361433997327, -15.7697199271455, 11.1373468568838,
	     -4.91933095295458, 1.39443093325178, -0.247689078940725,
	     0.0251861553415515, -0.00112130382664914}};

	double end[6] = {18.5527638190955, 13.7587939698492,
			 4.89447236180905,
			 11.4974874371859, 8.15075376884422,
			 3.17587939698492};
	int j;
	double c;
	switch((int)k[0]) {
	default: error("rho_ggw(): case (%i) not implemented.", (int)k[0]);
	    // c : identical numbers to those in SET_ABC_GGW  below
	case 1: j = 0; c = 1.694;     break;
	case 2: j = 1; c = 1.2442567; break;
	case 3: j = 2; c = 0.4375470; break;
	case 4: j = 3; c = 1.063;     break;
	case 5: j = 4; c = 0.7593544; break;
	case 6: j = 5; c = 0.2959132; break;
	}
	x = fabs(x);
	if (x <= c)
	    return(C[j][0]*x*x);
	else if (x <= 3*c)
	    return(C[j][1] +
		   x*(C[j][2] +
		      x*(C[j][3] +
			 x*(C[j][4] +
			    x*(C[j][5] +
			       x*(C[j][6] +
				  x*(C[j][7] +
				     x*(C[j][8] +
					x*(C[j][9])))))))));
	else if (x <= end[j])
	    return(C[j][10] +
		   x*(C[j][11] +
		      x*(C[j][12] +
			 x*(C[j][13] +
			    x*(C[j][14] +
			       x*(C[j][15] +
				  x*(C[j][16] +
				     x*(C[j][17] +
					x*(C[j][18]+
					   x*(C[j][19]))))))))));
	else return(1.);
    }
    else { // k[0] == 0; k[1:4] = (a, b, c, rho(Inf)) =  "general parameters"
	// calculate integral
	x = fabs(x);
	double a = 0., epsabs = R_pow(DBL_EPSILON, 0.25), result, abserr;
	int neval, ier, last, limit = 100, lenw = 4 * limit;
	int   *iwork =    (int *) R_alloc(limit, sizeof(int));
	double *work = (double *) R_alloc(lenw,  sizeof(double));
	Rdqags(psi_ggw_vec, (void *)k, &a, &x, &epsabs, &epsabs,
	       &result, &abserr, &neval, &ier,
	       &limit, &lenw, &last,
	       iwork, work);
	if (ier >= 1) {
	    error("Error while calling Rdqags(): ier = %i", ier);
	}
	return(result/k[4]);
    }
}

void psi_ggw_vec(double *x, int n, void *k)
{
    for (int i = 0; i<n; i++) x[i] = psi_ggw(x[i], k);
}

double psi_ggw(double x, const double k[])
{
    /*
     * Gauss Weight with constant center
     */
#define SET_ABC_GGW(NAME)					\
    /* set a,b,c */						\
	double a, b, c;						\
	switch((int)k[0]) {					\
	default: error(#NAME "_ggw: Case not implemented.");	\
	    /* user specified: */				\
	case 0: a = k[1];      b = k[2]; c = k[3]; break;	\
	    /* Set of predefined cases: */			\
	case 1: a = 0.648;     b = 1.;   c = 1.694; break;	\
	case 2: a = 0.4760508; b = 1.;   c = 1.2442567; break;	\
	case 3: a = 0.1674046; b = 1.;   c = 0.4375470; break;	\
	case 4: a = 1.387;     b = 1.5;  c = 1.063; break;	\
	case 5: a = 0.8372485; b = 1.5;  c = 0.7593544; break;	\
	case 6: a = 0.2036741; b = 1.5;  c = 0.2959132; break;	\
	}							\
	double ax = fabs(x);

    SET_ABC_GGW(psi);
    if (ax < c) return x;
    else {
	a = -R_pow(ax-c,b)/2/a;
	return (a < MIN_Exp) ? 0. : x*exp(a);
    }
}

double psip_ggw(double x, const double k[])
{
    /*
     * Gauss Weight with constant center
     */
    SET_ABC_GGW(psip);
    if (ax < c) return 1.;
    else {
	double ea;
	a *= 2.;
	ea = -R_pow(ax-c,b)/a;
	return (ea < MIN_Exp) ? 0. : exp(ea) * (1 - b/a*ax*R_pow(ax-c,b-1));
    }
}

double wgt_ggw(double x, const double k[])
{
    /*
     * Gauss Weight with constant center
     */
    SET_ABC_GGW(wgt);
    return (ax < c) ? 1. : exp(-R_pow(ax-c,b)/2/a);
}
#undef SET_ABC_GGW


//--- LQQ := Linear-Quadratic-Quadratic ("lqq") --------------------------------
//--- ---    was LGW := "lin psip" := piecewise linear psi'() ------------------

// k[0:2] == (b, c, s) :
// k[0]= b = bend adjustment
// k[1]= c = cutoff of central linear part
// k[2]= s : "slope of descending": 1 - s = min_x psi'(x)

// "lin psip" := piecewise linear psi'() :
double psip_lqq (double x, const double k[])
{
    double ax = fabs(x);
    if (ax <= k[1])
	return(1.);
    else {
	double k01 = k[0] + k[1];// = b+c
	if (/*k[1] < ax && */ ax <= k01)
	    return 1. - k[2]/k[0] * (ax - k[1]);
	else {
	    double
		s5 = 1. - k[2], // = (1-s)
 		a = (k[0] * k[2] - 2 * k01)/ s5;
	    if (/* k01 < ax && */ ax < k01 + a)
		return -s5*((ax - k01)/a -1.);
	    else
		return 0.;
	}
    }
}

// piecewise linear psi'()  ==> piecewise constant  psi''():
double psi2_lqq (double x, const double k[])
{
    // double sx = sign(x), ax = fabs(x); :
    double sx, ax;
    if (x < 0) { sx = -1; ax = -x; } else { sx = +1; ax = x; }

    // k[0:2] == (b, c, s) :
    if (ax <= k[1])
	return 0.;
    else {
	double k01 = k[0] + k[1];
	if (/*k[1] < ax && */ ax <= k01)
	    return sx * (- k[2]/k[0]);
	else {
	    double
		s5 = 1. - k[2], // = (1-s)
 		a = (k[0] * k[2] - 2 * k01)/ s5;
	    if (/* k01 < ax && */ ax < k01 + a)
		return sx * (- s5 / a);
	    else
		return 0.;
	}
    }
}

double psi_lqq (double x, const double k[])
{
    double ax = fabs(x);
    if (ax <= k[1])
	return(x);
    else {
	// k[0:2] == (b, c, s) :
	double k01 = k[0] + k[1];
	if (ax <= k01)
	    return((double) (x>0 ? 1 : (x<0 ? -1 : 0)) *
		   (ax - k[2] * pow(ax - k[1], 2.) / k[0] / 2.));
	else {
	    double
		s5 = k[2] - 1., // s - 1
		s6 = -2 * k01 + k[0] * k[2]; // numerator( -a ) ==> s6/s5 = -a
	    if (/* k01 < ax && */ ax < k01 - s6 / s5)
		return((double) (x>0 ? 1 : -1) *
		       (-s6/2. - pow(s5, 2.) / s6 * (pow(ax - k01, 2.) / 2. + s6 / s5 * (ax - k01))));
	    else
		return 0.;
	}
    }
}

double rho_lqq (double x, const double k[])
{
    double ax = fabs(x), k01 = k[0] + k[1];
    if (ax <= k[1])
	return((3. * k[2] - 3.) /
	       (k[2] * k[1] * (3. * k[1] + 2. * k[0]) +
		pow(k01, 2.)) * x * x);
    else if (/* k[1] < ax && */ ax <= k01) {
	double s0 = ax - k[1];
	return((6. * k[2] - 6.) /
	       (k[2] * k[1] * (3. * k[1] + 2. * k[0]) + pow(k01, 2.)) *
	       (x * x / 2. - k[2] / k[0] * pow(s0, 3.) / 6.));
    }
    else {
	double
	    s5 = k[2] - 1.,
	    s6 = -2 * k01 + k[0] * k[2];
	if (/* k01 < ax && */ ax < k01 - s6 / s5) {
	    double s7 = ax - k01, k01_2 = pow(k01, 2.);
	    return((6. * s5) /
		   (k[2] * k[1] * (3. * k[1] + 2. * k[0]) + k01_2) *
		   (k01_2 / 2. - k[2] * k[0] * k[0] / 6. -
		    s7/2. * (s6 + s7 * (s5 + s7 * s5 * s5 / 3. / s6))));
	}
	else
	    return 1.;
    }
}

double wgt_lqq (double x, const double k[])
{
    double ax = fabs(x);
    if (ax <= k[1])
	return(1.);
    else {
	double k01 = k[0] + k[1];
	if (ax <= k01) {
	    double s0 = ax - k[1];
	    return(1. - k[2] * s0 * s0 / (2 * ax * k[0]));
	}
	else {
	    double
		s5 = k[2] - 1.,
		s6 = -2 * k01 + k[0] * k[2];
	    if (ax < k01 - s6 / s5) {
		double s7 = ax - k01;
		return(-(s6/2. + s5 * s5 / s6 * s7 * (s7/2. + s6 / s5)) / ax);
	    }
	    else
		return(0.);
	}
    }
}

double rho_modOpt(double x, const double c[])
{
  /*
   * c[] contains 16 values:
   *   c[0]: tuning constant a
   *   c[1] & c[2]: endpoints of support interval calculated from a (x > 0)
   *   c[3]: rescaling parameter; 1.0 gives definition of optimal psi
   *   c[4]: Psi_opt(c[1], c)
   *   c[5]: rho_opt(Inf) when c[3] = 1
   *   c[6:10]: polynomial coefficients
   *   c[11:12]: endpoints of support (for x > 0)
   *   c[13]: rescaling parameter, 1.0 yields poly optimal psi
   *   c[14], c[15]: reference points (max value is c[15] - c[14])
   */
  
  x = fabs(x) / c[13];
  
  if(x <= c[11])
    return(0.0);
  
  if(x >= c[12])
    return(c[15]-c[14]);
  
  double x2 = x * x;
  
  double u3 = c[6] * x + c[7] * x2 / 2.0 + c[8] * x2 * x2 / 4.0 + c[9] * x2 * x2 * x2 / 6.0 + c[10] * x2 * x2 * x2 * x2 / 8.0;
  
  return(u3 - c[14]);  
}

double rho_modOptv0(double x, const double c[])
{
  /*
   * c[] contains 6 values:
   *   c[0]: tuning constant a
   *   c[1]: dnorm(1) / (dnorm(1) - a)
   *   c[2]: upper endpoint of optimal support interval
   *   c[3]: rescaling parameter; 1.0 gives definition of modified optimal psi
   *   c[4]: Psi_opt(1.0, c)
   *   c[5]: rho_modOpt(Inf) when c[3] = 1
   */
  
  x = fabs(x) / c[3];
  
  if(x < 1.0)
    return((0.5 * x * x) / c[5]);
  
  if(x > c[2])
    return(1.0);
  
  return((0.5 + c[1] * (Psi_opt(x, c) - c[4])) / c[5]);
}



double psi_modOpt(double x, const double c[])
{
 /*
  * c[] contains 6 values:
  *   c[0]: tuning constant a
  *   c[1]: dnorm(1) / (dnorm(1) - a)
  *   c[2]: upper endpoint of optimal support interval
  *   c[3]: rescaling parameter; 1.0 gives definition of modified optimal psi
  *   c[4]: Psi_opt(1.0, c)
  *   c[5]: rho_modOpt(Inf) when c[3] = 1
  */

  double absx = 0.0, y = 0.0;

  x = x / c[3];
  absx = fabs(x);

  if(absx <= 1.0)
    return(c[3] * x);

  if(absx >= c[2])
    return(0.0);

  y = c[3] * c[1] * (absx - c[0] / dnorm(absx, 0.0, 1.0, 0));
  return(x > 0.0 ? y : -1.0 * y);
}

double psip_modOpt(double x, const double c[])
{
 /*
  * c[] contains 6 values:
  *   c[0]: tuning constant a
  *   c[1]: dnorm(1) / (dnorm(1) - a)
  *   c[2]: upper endpoint of optimal support interval
  *   c[3]: rescaling parameter; 1.0 gives definition of modified optimal psi
  *   c[4]: Psi_opt(1.0, c)
  *   c[5]: rho_modOpt(Inf) when c[3] = 1
  */

  x = fabs(x) / c[3];

  if(x <= 1.0)
    return(1.0);

  if(x >= c[2])
    return(0.0);

  return(c[1] * (1.0 - c[0] * x / dnorm(x, 0.0, 1.0, 0)));
}

double wgt_modOpt(double x, const double c[])
{
 /*
  * c[] contains 6 values:
  *   c[0]: tuning constant a
  *   c[1]: dnorm(1) / (dnorm(1) - a)
  *   c[2]: upper endpoint of optimal support interval
  *   c[3]: rescaling parameter; 1.0 gives definition of modified optimal psi
  *   c[4]: Psi_opt(1.0, c)
  *   c[5]: rho_modOpt(Inf) when c[3] = 1
  */
  
  x = fabs(x) / c[3];

  if(x <= 1.0)
    return(1.0);

  if(x >= c[2])
    return(0.0);
  
  return(c[1] * (1.0 - c[0] / (x * dnorm(x, 0.0, 1.0, 0))));
}

/*============================================================================*/


/* this function finds the k-th place in the
 * vector a, in the process it permutes the
 * elements of a
 */
double kthplace(double *a, int n, int k)
{
    int jnc,j;
    int l,lr;
    double ax,w;
    k--;
    l=0;
    lr=n-1;
    while (l < lr) {
	ax=a[k];
	jnc=l;
	j=lr;
	while (jnc <= j) {
	    while (a[jnc] < ax) jnc++;
	    while (a[j] > ax) j--;
	    if (jnc <= j) {
		w=a[jnc];
		a[jnc]=a[j];
		a[j]=w;
		jnc++;
		j--;
	    }
	}
	if (j < k) l=jnc;
	if (k < jnc) lr=j;
    }
    return(a[k]);
}

/* This is from VR's bundle, MASS package  VR/MASS/src/lqs.c : */
/*
   Sampling k from 0:n-1 without replacement.
 */
static void sample_noreplace(int *x, int n, int k, int *ind_space)
{
    int i, j, nn=n;
#define II ind_space

    for (i = 0; i < n; i++) II[i] = i;
    for (i = 0; i < k; i++) {
	j = nn * unif_rand();
	x[i] = II[j];
	II[j] = II[--nn];
    }
#undef II
}

/* RWLS iterations starting from i_estimate,
 * ---- the workhorse of the "lmrob_MM" algorithm;
 * in itself,  ``just'' an M-estimator :
 */
Rboolean rwls(const double X[], const double y[], int n, int p,
	 double *estimate, double *i_estimate,
	 double *resid, double* loss,
	 double scale, double epsilon,
	 int *max_it, /* on Input:  maximal number of iterations;
			 on Output: number of iterations */
	 double *rho_c, const int ipsi, int trace_lev)
{
    int lwork = -1, one = 1, info = 1;
    double work0, *work, wtmp, *weights;
    double *wx, *wy, done = 1., dmone = -1.;
    double *beta0, d_beta = 0.;
    int j, k, iterations = 0;
    Rboolean converged = FALSE;

    wx     = (double *) R_alloc(n*p, sizeof(double));
    wy     = (double *) R_alloc(n,   sizeof(double));
    beta0  = (double *) R_alloc(p,   sizeof(double));

    INIT_WLS(wx, wy, n, p);

    COPY(i_estimate, beta0, p);
    /* calculate residuals */
    COPY(y, resid, n);
    F77_CALL(dgemv)("N", &n, &p, &dmone, X, &n, beta0, &one, &done, resid, &one FCONE);

    /* main loop */
    while(!converged &&	 ++iterations < *max_it) {
	R_CheckUserInterrupt();
        /* compute weights */
	get_weights_rhop(resid, scale, n, rho_c, ipsi, weights);
	/* solve weighted least squares problem */
	COPY(y, wy, n);
	FIT_WLS(X, wx, wy, n, p);
	COPY(wy, estimate, p);
	/* calculate residuals */
	COPY(y, resid, n);
	F77_CALL(dgemv)("N", &n, &p, &dmone, X, &n, estimate, &one, &done, resid, &one FCONE);
	if(trace_lev >= 3) {
	    /* get the residuals and loss for the new estimate */
	    *loss = sum_rho_sc(resid,scale,n,0,rho_c,ipsi);
	    Rprintf("  it %4d: L(b1) = %.12g ", iterations, *loss);
	}
	/* check for convergence */
	d_beta = norm1_diff(beta0,estimate, p);
	if(trace_lev >= 3) {
	    if(trace_lev >= 4) {
		Rprintf("\n  b1 = (");
		for(j=0; j < p; j++)
		    Rprintf("%s%.11g", (j > 0)? ", " : "", estimate[j]);
		Rprintf(");");
	    }
	    Rprintf(" ||b0 - b1||_1 = %g\n", d_beta);
	}
	converged = d_beta <= epsilon * fmax2(epsilon, norm1(estimate, p));
	COPY(estimate, beta0, p);
    } /* end while(!converged & iter <=...) */

    if (trace_lev < 3)
	/* get the residuals and loss for the new estimate */
	*loss = sum_rho_sc(resid,scale,n,0,rho_c,ipsi);

    if(trace_lev)
	Rprintf(" rwls() used %d it.; last ||b0 - b1||_1 = %g;%sconvergence\n",
		iterations, d_beta, (converged ? " " : " NON-"));

    *max_it = iterations;

    CLEANUP_WLS;

    return converged;
} /* rwls() */

/* sets the entries of a matrix to zero */
void zero_mat(double **a, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++)
	for(j=0; j < m; j++)
	    a[i][j] = 0.;
}

/*
 *
 * 2004 / 5 -- Matias Salibian-Barrera & Victor Yohai
 * Department of Statistics, University of British Columbia
 * matias@stat.ubc.ca
 * Department of Mathematics, University of Buenos Aires
 * vyohai@uolsinectis.com.ar
 *
 *
 * Reference: A fast algorithm for S-regression estimates,
 * 2005, Salibian-Barrera and Yohai.
 */

/* This function implements the "large n" strategy
 */
void fast_s_large_n(double *X, double *y,
		    int *nn, int *pp, int *nRes, int *max_it_scale,
		    int *ggroups, int *nn_group,
		    int *K, int *max_k, double rel_tol, double inv_tol, int *converged,
		    int *best_r, double *bb, double *rrhoc, int *iipsi,
		    double *bbeta, double *sscale,
		    int trace_lev, int mts, int ss)
{
/* *X  = the n x p  design matrix (incl. intercept if appropriate),
 *	 in column order as sent by R)
 * *y  = the ( n ) response vector
 * *nn =: n = the length of y
 * *pp =: p = the number of columns in X
 * *nRes  = number of re-sampling candidates to be used in each partition
 * *K     = number of refining steps for each candidate (typically 1 or 2)
 * *max_k = number of refining steps for each candidate (typically 1 or 2)
             [used to be hard coded to MAX_ITER_REFINE_S = 50 ]
 * *rel_tol= convergence tolerance for iterative refinement iterations
             [used to be hard coded to EPS = 1e-7 ]
 * *converged: will become 0(FALSE)  iff at least one of the best_r iterations
 *             did not converge (in max_k steps to rel_tol precision)
 * *ggroups = number of groups in which to split the
 *	      random subsample
 * *nn_group = size of each of the (*ggroups) groups
 *	       to use in the random subsample
 * *best_r = no. of best candidates to be iterated further ("refined")
 * *bb	   = right-hand side of S-equation (typically 1/2)
 * *rrhoc  = tuning constant for loss function
 *	     (this should be associated with *bb)
 * *iipsi  = indicator for type of psi function to be used
 * *bbeta  = final estimator
 * *sscale = associated scale estimator (or -1 when problem)
 */
    int i,j,k, k2, it_k, ij, freedsamp = 0, initwls = 0;
    int n = *nn, p = *pp, kk = *K, ipsi = *iipsi;
    int groups = *ggroups, n_group = *nn_group, sg = groups * n_group;
    double b = *bb, sc, best_sc, worst_sc;
    int pos_worst_scale;
    Rboolean conv;
    /* (Pointers to) Arrays - to be allocated */
    int *indices, *ind_space;
    double **best_betas, *best_scales;
    double *xsamp, *ysamp, *beta_ref;
    double **final_best_betas, *final_best_scales;

#define CALLOC_MAT(_M_, _n_, _d_)			\
    _M_ = (double **) Calloc(_n_, double *);		\
    for(int i=0; i < _n_; i++)				\
	_M_[i] = (double *) Calloc(_d_, double)

    beta_ref = (double *) Calloc(p, double);
    CALLOC_MAT(final_best_betas, *best_r, p);
    final_best_scales = (double *) Calloc(*best_r, double);
    k = *best_r * groups;
    best_scales = (double *) Calloc(k,	double );
    CALLOC_MAT(best_betas, k, p);
    indices =   (int *)    Calloc(sg,   int);
    ind_space = (int *)    Calloc(n,   int);
    xsamp =     (double *) Calloc(n_group*p, double);
    ysamp =     (double *) Calloc(n_group,   double);

    /* assume that n > 2000 */

    /*	set the seed */
    GetRNGstate();

    /* get a sample of k indices */
    sample_noreplace(indices, n, sg, ind_space);
    /* FIXME: define groups using nonsingular subsampling? */
    /*        would also need to allow observations to be part */
    /*        of multiple groups at the same time */
    Free(ind_space);
    /* FIXME: Also look at lqs_setup(),
     * -----  and  xr[.,.] "fortran-like" matrix can be used from there!*/

/* For each (of 'groups') group : get the *best_r best betas : */

#define X(_k_, _j_) X[_j_*n+_k_]
#define xsamp(_k_, _j_) xsamp[_j_*n_group+_k_]

    for(i=0; i < groups; i++) {
	/* populate matrix */
	for(j = 0; j < n_group; j++) {
	    ij = i*n_group + j;
	    for (k = 0; k < p; k++)
		xsamp(j, k) = X(indices[ij], k);
	    ysamp[j] = y[indices[ij]];
	}
	if (trace_lev)
	    Rprintf(" Subsampling to find candidate betas in group %d:\n", i);
	if(fast_s_with_memory(xsamp, ysamp,
			      &n_group, pp, nRes, max_it_scale, K, max_k, rel_tol, inv_tol,
			      trace_lev, best_r, bb, rrhoc,
			      iipsi, best_betas + i* *best_r,
			      best_scales+ i* *best_r, mts, ss)) {
	    *sscale = -1.; /* problem */
	    goto cleanup_and_return;
	}
    }

    Free(xsamp); Free(ysamp); freedsamp = 1;
#undef xsamp

/* now	iterate (refine) these "best_r * groups"
 * best betas in the (xsamp,ysamp) sample
 * with kk C-steps and keep only the "best_r" best ones
 */
    /* initialize new work matrices */
    double *wx, *wy, *res;
    res = (double *) R_alloc(n,   sizeof(double));
    wx =  (double *) R_alloc(n*p, sizeof(double)); // need only k here,
    wy =  (double *) R_alloc(n,   sizeof(double)); // but n in the last step
    xsamp =     (double *) Calloc(sg*p, double);
    ysamp =     (double *) Calloc(sg,   double);
    freedsamp = 0;

#define xsamp(_k_,_j_) xsamp[_j_*sg+_k_]

    for (ij = 0; ij < sg; ij++) {
	for (k = 0; k < p; k++)
	    xsamp(ij, k) = X(indices[ij],k);
	ysamp[ij] = y[indices[ij]];
    }

    int lwork = -1, one = 1, info = 1;
    double work0, *work, *weights;
    INIT_WLS(wx, wy, n, p); initwls = 1;

    conv = FALSE;
    pos_worst_scale = 0;
    for(i=0; i < *best_r; i++)
	final_best_scales[i] = INFI;
    worst_sc = INFI;
    /* set the matrix to zero */
    zero_mat(final_best_betas, *best_r, p);
    for(i=0; i < (*best_r * groups); i++) {
	if(trace_lev >= 3) {
	    Rprintf("  Sample[%3d]: before refine_(*, conv=FALSE):\n", i);
	    Rprintf("   beta_cand : "); disp_vec(best_betas[i],p);
	    Rprintf("   with scale %.15g\n", best_scales[i]);

	}
	refine_fast_s(xsamp, wx, ysamp, wy, weights, sg, p, res,
		      work, lwork, best_betas[i],
		      kk, &conv/* = FALSE*/, *max_k, rel_tol, trace_lev,
		      b, rrhoc, ipsi, best_scales[i], /* -> */ beta_ref, &sc);
	if(trace_lev >= 3) {
	    Rprintf("   after refine: beta_ref : "); disp_vec(beta_ref,p);
	    Rprintf("   with scale %.15g\n", sc);
	}
	if ( sum_rho_sc(res, worst_sc, sg, p, rrhoc, ipsi) < b ) {
	    /* scale will be better */
	    sc = find_scale(res, b, rrhoc, ipsi, sc, sg, p, *max_it_scale);
	    k2 = pos_worst_scale;
	    final_best_scales[ k2 ] = sc;
	    COPY(beta_ref, final_best_betas[k2], p);
	    pos_worst_scale = find_max(final_best_scales, *best_r);
	    worst_sc = final_best_scales[pos_worst_scale];
	}
    }

    Free(xsamp); Free(ysamp); freedsamp = 1;

/* now iterate the best "best_r"
 * betas in the whole sample until convergence (max_k, rel_tol)
 */
    best_sc = INFI; *converged = 1;  k = 0;
    if(trace_lev)
	Rprintf(" Now refine() to convergence for %d very best ones:\n",
		*best_r);

    for(i=0; i < *best_r; i++) {
	conv = TRUE;
	it_k = refine_fast_s(X, wx, y, wy, weights, n, p, res,
			     work, lwork, final_best_betas[i], kk,
			     &conv/* = TRUE */, *max_k, rel_tol, trace_lev,
			     b, rrhoc, ipsi, final_best_scales[i],
			     /* -> */ beta_ref, &sc);
	if(trace_lev)
	    Rprintf("  Best[%d]: %sconvergence (%d iter.)",
		    i, conv ? "" : "NON ", it_k);
	if(best_sc > sc) {
	    if(trace_lev)
		Rprintf(": -> improved scale to %.15g", sc);
	    best_sc = sc;
	    COPY(beta_ref, bbeta, p);
	}
	if (trace_lev) Rprintf("\n");
	if (!conv && *converged) *converged = 0;
	if (k < it_k) k = it_k;
    }
    *sscale = best_sc;
    *max_k = k;

/* Done. Now clean-up. */

  cleanup_and_return:
    PutRNGstate();

    Free(best_scales);
    k = *best_r * groups;
    for(i=0; i < k; i++) Free( best_betas[i] );
    Free(best_betas); Free(indices);
    for(i=0; i < *best_r; i++)
	Free(final_best_betas[i]);
    Free(final_best_betas);
    Free(final_best_scales);
    Free(beta_ref);

    if (freedsamp == 0) {
	Free(xsamp); Free(ysamp);
    }

    if (initwls) {
	CLEANUP_WLS;
    }

#undef X
#undef xsamp

} /* fast_s_large_n() */

int fast_s_with_memory(double *X, double *y,
		       int *nn, int *pp, int *nRes, int *max_it_scale,
		       int *K, int *max_k, double rel_tol, double inv_tol,
		       int trace_lev, int *best_r, double *bb, double *rrhoc,
		       int *iipsi, double **best_betas, double *best_scales,
		       int mts, int ss)
{
/*
 * Called from fast_s_large_n(), the adjustment for large "n",
 * same as fast_s, but it returns the best_r best betas,
 * and their associated scales.
 *
 * x an	 n x p design matrix (including intercept if appropriate)
 * y an	 n vector
 * *nn = n, *pp = p
 * *nRes   = number of re-sampling candidates to be taken
 * *K	   = number of refining steps for each candidate
 * *best_r = number of (refined) to be retained for full iteration
 * *bb	   = right-hand side of the S-equation (typically 1/2)
 * *rrhoc  = tuning constant for loss function
 *	     (this should be associated with *bb)
 * *iipsi  = indicator for type of loss function to be used
 * *best_betas	= returning the best ... coefficient vectors
 * *best_scales = returning their associated residual scales
 */

    int i,j,k;
    int n = *nn, p = *pp, nResample = *nRes;
    Rboolean conv = FALSE;
    double *beta_cand, *beta_ref, *res;
    int ipsi = *iipsi;
    double b = *bb, sc, worst_sc = INFI;
    double work0, *weights, *work, *wx, *wy;
    int lwork = -1, one = 1, info = 1;
    int pos_worst_scale, sing=0;

    SETUP_SUBSAMPLE(n, p, X, 1);
    INIT_WLS(X, y, n, p);

    res	=       (double *) Calloc(n,   double);
    wx =        (double *) Calloc(n*p, double);
    wy =        (double *) Calloc(n,   double);
    beta_cand = (double *) Calloc(p,   double);
    beta_ref  = (double *) Calloc(p,   double);

    for(i=0; i < *best_r; i++)
	best_scales[i] = INFI;
    pos_worst_scale = 0;

/* resampling approximation  */

    for(i=0; i < nResample; i++) {
	R_CheckUserInterrupt();
	/* find a candidate */
	sing = subsample(Xe, y, n, p, beta_cand, ind_space, idc, idr, lu, v, pivot,
			 Dr, Dc, rowequ, colequ, 1, mts, ss, inv_tol, 1);
	if (sing) {
	    for (k=0; k< *best_r; k++) best_scales[i] = -1.;
	    goto cleanup_and_return;
	}
	/* FIXME: is_ok ?? */

	/* improve the re-sampling candidate */

	/* conv = FALSE : do *K refining steps */
	refine_fast_s(X, wx, y, wy, weights, n, p, res,
		      work, lwork, beta_cand, *K, &conv/* = FALSE*/,
		      *max_k, rel_tol, trace_lev, b, rrhoc, ipsi, -1.,
		      /* -> */ beta_ref, &sc);

	/* FIXME: if sc ~ 0 ---> return beta_cand and be done */

	if ( sum_rho_sc(res, worst_sc, n, p, rrhoc, ipsi) < b )	{
	    /* scale will be better */
	    sc = find_scale(res, b, rrhoc, ipsi, sc, n, p, *max_it_scale);
	    k = pos_worst_scale;
	    best_scales[ k ] = sc;
	    for(j=0; j < p; j++)
		best_betas[k][j] = beta_ref[j];
	    pos_worst_scale = find_max(best_scales, *best_r);
	    worst_sc = best_scales[pos_worst_scale];
	    if (trace_lev >= 2) {
	      Rprintf("  Sample[%3d]: found new candidate with scale %.7g\n",
		      i, sc);
	      Rprintf("               worst scale is now %.7g\n", worst_sc);
	    }
	}
    } /* for(i ) */

  cleanup_and_return:

    CLEANUP_SUBSAMPLE;
    CLEANUP_WLS;

    Free(res); Free(wx); Free(wy);
    Free(beta_cand); Free(beta_ref);

    return sing;
} /* fast_s_with_memory() */

void fast_s(double *X, double *y,
	    int *nn, int *pp, int *nRes, int *max_it_scale,
	    int *K, int *max_k, double rel_tol, double inv_tol, int *converged,
	    int *best_r, double *bb, double *rrhoc, int *iipsi,
	    double *bbeta, double *sscale, int trace_lev, int mts, int ss)
{
/* *X  = the n x p  design matrix (incl. intercept if appropriate),
 *	 in column order as sent by R)
 * *y  = the ( n ) response vector
 * *nn =: n = the length of y
 * *pp =: p = the number of columns in X
 * *nRes   = number of re-sampling candidates to be taken
 * *K	   = number of refining steps for each candidate
 * *best_r = number of (refined) to be retained for full iteration
 * *converged: will become FALSE  iff at least one of the best_r iterations
 *	       did not converge (in max_k steps to rel_tol precision)
 * *bb	   = right-hand side of the S-equation (typically 1/2)
 * *rrhoc  = tuning constant for loss function
 *	     (this should be associated with *bb)
 * *iipsi  = indicator for type of loss function to be used
 * *bbeta  = final estimator
 * *sscale = associated scale estimator (or -1 when problem)
 */
    int i,k, it_k;
    int n = *nn, p = *pp, nResample = *nRes, ipsi = *iipsi;
    Rboolean conv;
    double b = *bb;
    double sc, best_sc, worst_sc, aux;
    int pos_worst_scale;
    int lwork = -1, one = 1, info = 1;
    double work0, *work, *weights, sing=0;

    /* Rprintf("fast_s %d\n", ipsi); */

    /* (Pointers to) Arrays - to be allocated */
    double *wx, *wy, *beta_cand, *beta_ref, *res;
    double **best_betas, *best_scales;

    SETUP_SUBSAMPLE(n, p, X, 0);

    res	   = (double *) R_alloc(n, sizeof(double));
    wx     = (double *) R_alloc(n*p, sizeof(double));
    wy     = (double *) R_alloc(n,   sizeof(double));

    best_betas = (double **) Calloc(*best_r, double *);
    best_scales = (double *) Calloc(*best_r, double);
    for(i=0; i < *best_r; i++) {
	best_betas[i] = (double*) Calloc(p, double);
	best_scales[i] = INFI;
    }
    beta_cand = (double *) Calloc(p, double);
    beta_ref  = (double *) Calloc(p, double);

    INIT_WLS(wx, wy, n, p);

    /* disp_mat(x, n, p); */

    pos_worst_scale = 0; conv = FALSE; worst_sc = INFI;

    /* srand((long)*seed_rand); */
    GetRNGstate();

/* resampling approximation  */

    if (trace_lev)
	Rprintf(" Subsampling to find candidate betas:\n", i);

    for(i=0; i < nResample; i++) {

	R_CheckUserInterrupt();
	/* find a candidate */
	sing = subsample(Xe, y, n, p, beta_cand, ind_space, idc, idr, lu, v, pivot,
			 Dr, Dc, rowequ, colequ, 1, mts, ss, inv_tol, 1);
	if (sing) {
	    *sscale = -1.;
	    goto cleanup_and_return;
	}
	if (trace_lev >= 5) {
	    Rprintf("  Sample[%3d]: idc = ", i); disp_veci(idc, p);
	}

	/* disp_vec(beta_cand,p); */

	/* improve the re-sampling candidate */

	/* conv = FALSE : do *k refining steps */
	refine_fast_s(X, wx, y, wy, weights, n, p, res,
		      work, lwork, beta_cand, *K, &conv/* = FALSE*/,
		      *max_k, rel_tol, trace_lev, b, rrhoc, ipsi, -1.,
		      /* -> */ beta_ref, &sc);
	if(trace_lev >= 3) {
	    double del = norm_diff(beta_cand, beta_ref, p);
	    Rprintf("  Sample[%3d]: after refine_(*, conv=FALSE):\n", i);
	    Rprintf("   beta_ref : "); disp_vec(beta_ref,p);
	    Rprintf("   with ||beta_ref - beta_cand|| = %.12g, --> sc = %.15g\n",
		    del, sc);
	}
	if(fabs(sc) == 0.) { /* exact zero set by refine_*() */
	    if(trace_lev >= 1)
		Rprintf(" Too many exact zeroes -> leaving refinement!\n");
	    *sscale = sc;
	    COPY(beta_cand, bbeta, p);
	    goto cleanup_and_return;
	}
	if ( sum_rho_sc(res, worst_sc, n, p, rrhoc, ipsi) < b )	{
	    /* scale will be better */
	    sc = find_scale(res, b, rrhoc, ipsi, sc, n, p, *max_it_scale);
	    k = pos_worst_scale;
	    best_scales[ k ] = sc;
	    COPY(beta_ref, best_betas[k], p);
	    pos_worst_scale = find_max(best_scales, *best_r);
	    worst_sc = best_scales[pos_worst_scale];
	    if (trace_lev >= 2) {
	      Rprintf("  Sample[%3d]: found new candidate with scale %.7g\n", i, sc);
	      Rprintf("               worst scale is now %.7g\n", worst_sc);
	    }
	}

    } /* for(i ) */

/* now look for the very best */
    if(trace_lev)
	Rprintf(" Now refine() to convergence for %d very best ones:\n",
		*best_r);

    best_sc = INFI; *converged = 1;  k = 0;
    for(i=0; i < *best_r; i++) {
	conv = TRUE;
	if(trace_lev >= 4) Rprintf("  i=%d:\n", i);
	it_k = refine_fast_s(X, wx, y, wy, weights, n, p, res, work, lwork,
			     best_betas[i], *K,  &conv /* = TRUE */, *max_k,
			     rel_tol, trace_lev, b, rrhoc, ipsi,
			     best_scales[i], /* -> */ beta_ref, &aux);
	if(trace_lev)
	    Rprintf("  Best[%d]: %sconvergence (%d iter.)",
		    i, (conv) ? "" : "NON ", it_k);
	if(aux < best_sc) {
	    if(trace_lev)
		Rprintf(": -> improved scale to %.15g", aux);
	    best_sc = aux;
	    COPY(beta_ref, bbeta, p);
	}
	if(trace_lev) Rprintf("\n");
	if (!conv && *converged) *converged = 0;
	if (k < it_k) k = it_k;
    }
    *sscale = best_sc;
    *max_k = k;

  cleanup_and_return:

    PutRNGstate();

    CLEANUP_SUBSAMPLE;
    CLEANUP_WLS;

    Free(best_scales);
    Free(beta_cand);
    Free(beta_ref);
    for(i=0; i < *best_r; i++)
	Free(best_betas[i]);
    Free(best_betas);

    return;
} /* fast_s() */

int refine_fast_s(const double X[], double *wx, const double y[], double *wy,
		  double *weights,
		  int n, int p, double *res,
		  double *work, int lwork, double *beta_cand,
		  int kk, Rboolean *conv, int max_k, double rel_tol,
		  int trace_lev,
		  double b, double *rrhoc, int ipsi, double initial_scale,
		  double *beta_ref, double *scale)
{
/*
 * X	   = matrix (n x p) of explanatory variables
 * y	   = vector ( n )   of responses
 * weights = robustness weights wt[] * y[]	(of length n)
 * res	   = residuals	y[] - x[,] * beta	(of length n)
 * conv: FALSE means do kk refining steps      (and conv stays FALSE)
 *	 TRUE  means refine until convergence(rel_tol, max_k)
 *             and in this case, 'conv' *returns* TRUE if refinements converged
 * beta_cand= candidate beta[] (of length p)	Input *and* Output
 * is	    = initial scale			input

 * beta_ref = resulting beta[] (of length p)	Output
 * scale    = final scale			Output

 * for FIT_WLS, DGELS:
 * wx       = matrix (n x p)
 * wy       = vector of length n
 * work     = vector of length lwork
 * lwork    = length of vector work
 */

    int i,j,k, zeroes=0, one = 1, info = 1;
    Rboolean converged = FALSE;/* Wall */
    double s0, done = 1., dmone = -1., wtmp;

    if (trace_lev >= 4) {
	Rprintf("   beta_cand before refinement : "); disp_vec(beta_cand,p);
    }

    /* calculate residuals */
    COPY(y, res, n);
    F77_CALL(dgemv)("N", &n, &p, &dmone, X, &n, beta_cand, &one, &done, res, &one FCONE);
    for(j=0; j < n; j++) {
	if( fabs(res[j]) < EPS_SCALE )
	    zeroes++;
    }
/* if "perfect fit", return it with a 0 assoc. scale */
    if( zeroes > (((double)n + (double)p)/2.) ) /* <<- FIXME: depends on 'b' ! */
	{
	    COPY(beta_cand, beta_ref, p);
	    *scale = 0.;
	    return 0;
	}

    if( initial_scale < 0. )
	initial_scale = MAD(res, n, 0., wy, weights);// wy and weights used as work arrays
    s0 = initial_scale;
    if(*conv)
	kk = max_k;
    for(i=0; i < kk; i++) {
	/* one step for the scale */
	s0 = s0 * sqrt( sum_rho_sc(res, s0, n, p, rrhoc, ipsi) / b );
	/* compute weights for IRWLS */
	get_weights_rhop(res, s0, n, rrhoc, ipsi, weights);
        /* solve weighted least squares problem */
	COPY(y, wy, n);
	FIT_WLS(X, wx, wy, n, p);
	COPY(wy, beta_ref, p);
	if(*conv) { /* check for convergence */
	    double del = norm_diff(beta_cand, beta_ref, p);
	    double nrmB= norm(beta_cand, p);
	    if(trace_lev >= 4)
		Rprintf("   it %4d, ||b[i]||= %.12g, ||b[i] - b[i-1]|| = %.15g\n",
			i, nrmB, del);
	    converged = (del <= rel_tol * fmax2(rel_tol, nrmB));
	    if(converged)
		break;
	}
	/* calculate residuals */
	COPY(y, res, n);
	F77_CALL(dgemv)("N", &n, &p, &dmone, X, &n, beta_ref, &one, &done, res, &one FCONE);
	COPY(beta_ref, beta_cand, p);
    } /* for(i = 0; i < kk ) */

    if(*conv) {
	/* was "if(0)",	 since default lead to 'NOT converged' */
	if(!converged) {
	    *conv = FALSE;
	    warning("S refinements did not converge (to refine.tol=%g) in %d (= k.max) steps",
		    rel_tol, i);
	}
    }
    *scale = s0;
    return i; /* number of refinement steps */
} /* refine_fast_s() */

/* Subsampling part for M-S algorithm                    */
/* Recreates RLFRSTML function found in src/lmrobml.f    */
/* of the robust package                                 */
void m_s_subsample(double *X1, double *y, int n, int p1, int p2,
		   int nResample, int max_it_scale, double rel_tol, double inv_tol, double *bb,
		   double *rrhoc, int ipsi, double *sscale, int trace_lev,
		   double *b1, double *b2, double *t1, double *t2,
		   double *y_tilde, double *res, double *x1, double *x2,
		   int *NIT, int *K, int *KODE, double *SIGMA, double *BET0,
		   double *SC1, double *SC2, double *SC3, double *SC4, int mts, Rboolean ss)
{
    int i, one = 1, p = p1 + p2, info;
    double b = *bb, sc = INFI, done = 1., dmone = -1.;
    *sscale = INFI;

    if (trace_lev >= 2)
	Rprintf(" Starting subsampling procedure.. ");

    SETUP_SUBSAMPLE(n, p2, x2, 0);

    /*	set the seed */
    GetRNGstate();

    if (trace_lev >= 2) Rprintf(" [setup Ok]\n");

    for(i=0; i < nResample; i++) {
	R_CheckUserInterrupt();
	/* STEP 1: Draw a subsample of size p2 from (X2, y) */
	Rboolean sing = subsample(Xe, y, n, p2, t2, ind_space, idc, idr, lu, v, pivot,
				  Dr, Dc, rowequ, colequ, /* sample= */ TRUE, mts,
				  ss, inv_tol, /*solve = */ TRUE);
	if (sing) {
	    *sscale = -1.;
	    goto cleanup_and_return;
	}
	/* calculate partial residuals */
	COPY(y, y_tilde, n);
        F77_CALL(dgemv)("N", &n, &p2, &dmone, x2, &n, t2, &one, &done, y_tilde, &one FCONE);
	/* STEP 3: Obtain L1-estimate of b1 */
	COPY(X1, x1, n*p1);
	F77_CALL(rllarsbi)(x1, y_tilde, &n, &p1, &n, &n, &rel_tol,
			   NIT, K, KODE, SIGMA, t1, res, SC1, SC2,
			   SC3, SC4, BET0);
	if (*KODE > 1) {
	    REprintf("m_s_subsample(): Problem in RLLARSBI (RILARS). KODE=%d. Exiting.\n",
		     *KODE);
	    *sscale = -1.;
	    goto cleanup_and_return;
	}
	/* STEP 4: Check if candidate looks promising */
	if (sum_rho_sc(res, *sscale, n, p, rrhoc, ipsi) < b) {
	    /* scale will be better */
	    /* STEP 5: Solve for sc */
	    sc = find_scale(res, b, rrhoc, ipsi, sc, n, p, max_it_scale);
	    if(trace_lev >= 2)
		Rprintf("  Sample[%3d]: new candidate with sc = %10.5g\n",i,sc);
	    /* STEP 6: Update best fit */
	    *sscale = sc;
	    COPY(t1, b1, p1);
	    COPY(t2, b2, p2);
	    if (sc < EPS_SCALE) {
		REprintf("\nScale too small\n",
			 "Aborting m_s_subsample()\n\n");
		*sscale = -1.;
		goto cleanup_and_return;
	    }
	}
    } /* for(i ) */

    /* STEP 7: Clean up and return */
    if (trace_lev >= 1) {
	Rprintf(" Finished M-S subsampling with scale = %.5f\n",*sscale);
#define maybe_SHOW_b1_b2			\
	if (trace_lev >= 3) {			\
	     Rprintf("  b1: "); disp_vec(b1,p1);\
	     Rprintf("  b2: "); disp_vec(b2,p2);\
	}
	maybe_SHOW_b1_b2;
    }

  cleanup_and_return:
    CLEANUP_SUBSAMPLE;
    PutRNGstate();
} /* m_s_subsample() */

/* Descent step for M-S algorithm
 * Return value: convergence; note that convergence is *not* guaranteed
 */
Rboolean m_s_descent(double *X1, double *X2, double *y,
		 int n, int p1, int p2, int K_m_s, int max_k, int max_it_scale,
		 double rel_tol, double *bb, double *rrhoc,  int ipsi,
		 double *sscale, int trace_lev,
		 double *b1, double *b2, double *t1, double *t2,
		 double *y_tilde, double *res, double *res2, double *x1, double *x2,
		 int *NIT, int *K, int *KODE, double *SIGMA,  double *BET0,
		 double *SC1, double *SC2, double *SC3, double *SC4)
{
    int j, k, nnoimpr = 0, nref = 0;
    int p = p1 + p2;
    Rboolean converged = FALSE;
    double b = *bb;
    double sc = *sscale, done = 1., dmone = -1.;
    int lwork = -1, one = 1, info = 1;
    double work0, *work, wtmp, *weights;

    COPY(b1, t1, p1);
    COPY(b2, t2, p2);
    COPY(res, res2, n);

    if (trace_lev >= 2)
	Rprintf(" Starting descent procedure...\n");

    INIT_WLS(x2, y, n, p2);

    if (trace_lev >= 3) {
	Rprintf("  Scale: %.5f\n", *sscale);
	if (trace_lev >= 5) {
	    Rprintf("   res2: "); disp_vec(res2,n);
	}
    }

    /* Do descent steps until there is no improvement for   */
    /* K_m_s steps or we are converged                      */
    /* (convergence is not guaranteed)                      */
    while ( (nref++ < max_k) & (!converged) & (nnoimpr < K_m_s) ) {
	R_CheckUserInterrupt();
	/* STEP 1: update b2 (save it to t2) */
	/* y_tilde = y - x1 %*% t1 */
	COPY(y, y_tilde, n);
	COPY(X1, x1, n*p1);
	F77_CALL(dgemv)("N", &n, &p1, &dmone, x1, &n, t1, &one, &done, y_tilde, &one FCONE);
	/* compute weights */
	get_weights_rhop(res2, sc, n, rrhoc, ipsi, weights);
	/* solve weighted least squares problem */
	FIT_WLS(X2, x2, y_tilde, n, p2);
	COPY(y_tilde, t2, p2);
        /* get (intermediate) residuals */
	COPY(y, res2, n);
	F77_CALL(dgemv)("N", &n, &p2, &dmone, X2, &n, t2, &one, &done, res2, &one FCONE);
	/* STEP 2: Obtain L1-estimate of b1 */
	COPY(res2, y_tilde, n);
	F77_CALL(rllarsbi)(x1, y_tilde, &n, &p1, &n, &n, &rel_tol,
			   NIT, K, KODE, SIGMA, t1, res2,
			   SC1, SC2, SC3, SC4, BET0);
	if (*KODE > 1) {
	    CLEANUP_WLS;
	    error("m_s_descent(): Problem in RLLARSBI (RILARS). KODE=%d. Exiting.",
		  *KODE);
	}
	/* STEP 3: Compute the scale estimate */
	sc = find_scale(res2, b, rrhoc, ipsi, sc, n, p, max_it_scale);
	/* STEP 4: Check for convergence */
	/* FIXME: check convergence using scale ? */
	double del = sqrt(norm_diff2(b1, t1, p1) + norm_diff2(b2, t2, p2));
	double nrmB = sqrt(norm2(t1, p1) + norm2(t2, p2));
	converged = (del < rel_tol * fmax2(rel_tol, nrmB));
	if (trace_lev >= 3) {
	    if(converged) Rprintf(" -->> converged\n");
	    if (trace_lev >= 4) {
		Rprintf("   Ref.step %3d: #{no-improvements}=%3d; (del,dB)=(%12.7g,%12.7g)\n",
			nref, nnoimpr, del, rel_tol * fmax2(rel_tol, nrmB));
		if (trace_lev >= 5) {
		    Rprintf("    weights: "); disp_vec(weights,n);
		    Rprintf("    t2: "); disp_vec(t2,p2);
		    Rprintf("    t1: "); disp_vec(t1,p1);
		    Rprintf("    res2: "); disp_vec(res2,n);
		}
	    }
	}
	/* STEP 5: Update best fit */
	if (sc < *sscale) {
	    COPY(t1, b1, p1);
	    COPY(t2, b2, p2);
	    COPY(res2, res, n);
	    *sscale = sc;
	    if (trace_lev >= 2)
		Rprintf("  Refinement step %3d: better fit, scale: %10.5g\n",
			nref, sc);
	    nnoimpr = 0;
	} else {
	    if (trace_lev >= 3)
		Rprintf("  Refinement step %3d: no improvement, scale: %10.5g\n",
			nref, sc);
	    nnoimpr++;
	}
    } // while(.)

    if ( (!converged) & (nref == max_k) )
	warning(" M-S estimate: maximum number of refinement steps reached.");

    if (trace_lev >= 1) {
	Rprintf(" Descent procedure: %sconverged (best scale: %.5g, last step: %.5g)\n",
		converged ? "" : "not ", *sscale, sc);
	if (nnoimpr == K_m_s)
	    Rprintf("  The procedure stopped after %d steps because there was no improvement in the last %d steps.\n  To increase this number, use the control parameter 'k.m_s'.\n", nref, nnoimpr);
	else if (trace_lev >= 2)
	    Rprintf("  No improvements in %d out of %d steps.\n", nnoimpr, nref);
	maybe_SHOW_b1_b2;
    }

    CLEANUP_WLS;

    return converged;
} /* m_s_descent() */

/* draw a subsample of observations and calculate a candidate           *
 * starting value for S estimates                                       *
 * uses a custom LU decomposition, which acts on the transposed design  *
 * matrix. In case of a singular subsample, the subsample is modified   *
 * until it is non-singular (for ss == 1).                              *
 *                                                                      *
 * Parts of the algorithm are based on the Gaxpy version of the LU      *
 * decomposition with partial pivoting by                               *
 * Golub G. H., Van Loan C. F. (1996) - MATRIX Computations             */
Rboolean subsample(const double x[], const double y[], int n, int m,
		   double *beta, int *ind_space, int *idc, int *idr,
		   double *lu, double *v, int *pivot,
		   double *Dr, double *Dc, int rowequ, int colequ,
		   Rboolean sample, int mts, Rboolean ss, double tol_inv, Rboolean solve)
{
    /* x:         design matrix (n x m)
       y:         response vector
       n:         length of y, nrow of x
       m:         ncol of x  ( == p )
       beta:      [out] candidate parameters (length m)
       ind_space: (required in sample_noreplace, length n)
                  holds the index permutation
       idc:       (required in sample_noreplace, !! length n !!)
		  [out] index of observations used in subsample
       idr:       work array of length m
       lu:        [out] LU decomposition of subsample of xt (m x m)
                  Note: U has is not rescaled by 1 / *cf, as it should,
		        this is done R_subsample().
       v:         work array of length m
       pivot:     [out] pivoting table of LU decomposition (length m-1)
       Dr:        row equilibration (as calculated in SETUP_EQUILIBRATION)
       Dc:        column equilibration
       rowequ:    whether rows were equilibrated
       coleq:     whether cols were equilibrated
       sample:    whether to sample or not
       mts:       the number of singular samples allowed before
                  giving up (Max Try Samples)
       ss:        type of subsampling to be used:
                  0: simple subsampling
                  1: nonsingular subsampling
       tol_inv:   tolerance for declaring a matrix singular
       solve:     solve the least squares problem on the subsample?
                  (0: no, 1: yes)

       return condition:
             0: success
             1: singular (matrix xt does not contain a m dim. full rank
                          submatrix)
             2: too many singular resamples (simple subsampling case)    */
    int j, k, l, one = 1, mu = 0, tmpi, i = 0, attempt = 0;
    double tmpd;
    Rboolean sing;

#define xt(_k_, _j_) x[idr[_k_]*n+idc[_j_]]
#define U(_k_, _j_) lu[_j_*m+_k_]
#define u(_k_, _j_) lu + (_j_*m+_k_)
#define L(_k_, _j_) lu[_j_*m+_k_]
#define l(_k_, _j_) lu + (_j_*m+_k_)

Start:
    /* STEP 1: Calculate permutation of 1:n */
    if (sample) {
	sample_noreplace(ind_space, n, n, idc);
    } else for(k=0;k<n;k++) ind_space[k] = k;
    for(k=0;k<m;k++) idr[k] = k;

    /* STEP 2: Calculate LU decomposition of the first m cols of xt     *
     *         using the order in ind_space                             */
    for(j = 0; j < m; j++) {
	sing=TRUE;
	do {
	    if (i+j == n) {
		warning("subsample(): could not find non-singular subsample.");
		return(1);
	    }
	    idc[j] = ind_space[i+j];
	    if (j == 0) {
		for(k=j;k<m;k++) v[k] = xt(k, j);
	    } else {
		for(k=0;k<j;k++) U(k,j) = xt(k, j);
		/* z = solve(lu[0:(j-1), 0:(j-1)], xt[0:(j-1), j]) */
		F77_CALL(dtrsv)("L", "N", "U", &j, lu, &m, u(0, j), &one FCONE FCONE FCONE);
		/* Rprintf("Step %d: z = ", j);  */
		/* for(i=0; i < j; i++) Rprintf("%lf ",U(i, j)); */
		/* Rprintf("\n"); */
		/* v[j:(m-1)] = xt[j:(m-1), j] - L[j:(m-1), 0:(j-1)] %*% z */
		for(k=j;k<m;k++) {
		    v[k] = xt(k, j);
		    for(l=0;l<j;l++) v[k] -= L(k, l) * U(l, j);
		}
		/* Rprintf("v = "); disp_vec(v, m); */
	    }
	    if (j < m-1) {
		/* find pivot */
		tmpd=fabs(v[j]); mu = j;
		for(k=j+1;k<m;k++) if (tmpd < fabs(v[k])) { mu = k; tmpd = fabs(v[k]); }
                /* debug possumDiv example, see tests/subsample.R */
		/* if (j == 36) { */
		/*     Rprintf("Step %d: ", j+1); */
		/*     for(k=j;k<m;k++) Rprintf("%lf ", fabs(v[k])); */
		/*     Rprintf("\n %d %lf\n", mu+1, v[mu]); */
		/*     Rprintf("47 > 51: %x\n", fabs(v[46]) > fabs(v[50])); */
		/*     Rprintf("47 < 51: %x\n", fabs(v[46]) < fabs(v[50])); */
		/* } */
		/* continue only if pivot is large enough */
		if (tmpd >= tol_inv) {
		    pivot[j] = mu;
		    tmpd = v[j]; v[j] = v[mu]; v[mu] = tmpd;
		    tmpi = idr[j]; idr[j] = idr[mu]; idr[mu] = tmpi;
		    for(k=j+1;k<m;k++) L(k, j) = v[k] / v[j];
		    if (j > 0) {
			for(k=0;k<j;k++) {
			    tmpd = L(j, k); L(j, k) = L(mu, k); L(mu, k) = tmpd;
			}
		    }
		}
	    }
	    if (fabs(v[j]) < tol_inv) {
		if (ss == 0) {
		    attempt++;
		    if (attempt >= mts) {
			warning("Too many singular resamples. Aborting subsample().\n See parameter 'subsampling; in help of lmrob.config().");
			return(2);
		    }
		    goto Start;
		}
		/* drop observation and try next one */
		i++;
	    } else {
		sing = FALSE;
		U(j, j) = v[j];
	    }
	} while(sing);
    } /* end for loop */

    /* Rprintf("lu:"); disp_vec(lu, m*m); */
    /* Rprintf("pivot:"); disp_veci(pivot, m-1); */
    /* Rprintf("idc:"); disp_veci(idc, m); */

    /* STEP 3: Solve for candidate parameters if requested */
    if (solve == 0) {
      for(k=0;k<m;k++) beta[k] = NA_REAL;
    } else {
      for(k=0;k<m;k++) beta[k] = y[idc[k]];
      /* scale y ( = beta ) */
      if (rowequ) for(k=0;k<m;k++) beta[k] *= Dr[idc[k]];
      /* solve U\tr L\tr \beta = y[subsample] */
      F77_CALL(dtrsv)("U", "T", "N", &m, lu, &m, beta, &one FCONE FCONE FCONE);
      F77_CALL(dtrsv)("L", "T", "U", &m, lu, &m, beta, &one FCONE FCONE FCONE);
      /* scale the solution */
      if (colequ) for(k=0;k<m;k++) beta[k] *= Dc[idr[k]];
      /* undo pivoting */
      for(k=m-2;k>=0;k--) {
    	tmpd = beta[k]; beta[k] = beta[pivot[k]]; beta[pivot[k]] = tmpd;
      }
    }

    return(0);

#undef Xt
#undef U
#undef u
#undef L
#undef l
}

void get_weights_rhop(const double r[], double s, int n, const double rrhoc[], int ipsi,
		      double *w)
{
    for(int i=0; i < n; i++)
	w[i] = wgt(r[i] / s, rrhoc, ipsi);
}

double find_scale(double *r, double b, double *rrhoc, int ipsi,
		  double initial_scale, int n, int p, int max_iter)
{
    double scale = initial_scale;
    for(int it = 0; it < max_iter; it++) {
	scale = initial_scale *
	    sqrt( sum_rho_sc(r, initial_scale, n, p, rrhoc, ipsi) / b ) ;
	if(fabs(scale - initial_scale) <= EPS_SCALE*initial_scale) // converged:
	    return(scale);
	initial_scale = scale;
    }
    warning("find_scale() did not converge in '%s' (= %d) iterations",
	    "maxit.scale", /* <- name from lmrob.control() */ max_iter);
    return(scale);
}

int find_max(double *a, int n)
{
    if(n==1)
	return(0);
    else {
	int i, k = 0;
	double tt = a[0];
	for(i=1; i < n; i++)
	    if(tt < a[i]) {
		tt = a[i];
		k = i;
	    }
	return(k);
    }
}

double sum_rho_sc(const double r[], double scale, int n, int p, const double c[], int ipsi)
{
    double s = 0;
    for(int i=0; i < n; i++)
	s += rho(r[i]/scale, c, ipsi);
    return(s / ((double) n - p));
}

/* ||x||_2^2 */
double norm2(double *x, int n)
{
    double s = 0.;
    int one = 1;
    s = F77_CALL(dnrm2)(&n, x, &one);
    return( s*s );
}

/* ||x||_2 */
double norm(double *x, int n)
{
    int one = 1;
    return(F77_CALL(dnrm2)(&n, x, &one));
}

double norm1(double *x, int n)
{
    int one = 1;
    return(F77_CALL(dasum)(&n, x, &one));
}

/* ||x-y||_2^2 */
double norm_diff2(double *x, double *y, int n)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
	s += (x[i]-y[i])*(x[i]-y[i]);
    return( s );
}

/* ||x-y||_2 */
double norm_diff(double *x, double *y, int n)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
	s += (x[i]-y[i])*(x[i]-y[i]);
    return( sqrt(s) );
}


/* ||x-y||_1 */
double norm1_diff(double *x, double *y, int n)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
	s += fabs(x[i]-y[i]);
    return(s);
}


double MAD(double *a, int n, double center, double *b,
	   double *tmp)
{
/* if center == 0 then do not center */
    int i;
/*     if( fabs(center) > 0.) { */
	for(i=0; i < n; i++)
	    b[i] = a[i] - center;
/*     } */
    return( median_abs(b,n,tmp) * 1.4826 );
}

double median(double *x, int n, double *aux)
{
    double t;
    int i;
    for(i=0; i < n; i++) aux[i]=x[i];
    if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2.0 ;
    else t = kthplace(aux,n, n/2+1 ) ;
    return(t);
}

double median_abs(double *x, int n, double *aux)
{
    double t;
    int i;
    for(i=0; i < n; i++) aux[i]=fabs(x[i]);
    if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2.0 ;
    else	t = kthplace(aux,n, n/2+1 ) ;
    return(t);
}

void disp_vec(double *a, int n)
{
    int i;
    for(i=0; i < n; i++) Rprintf("%lf ",a[i]);
    Rprintf("\n");
}

void disp_veci(int *a, int n)
{
    int i;
    for(i=0; i < n; i++) Rprintf("%d ",a[i]);
    Rprintf("\n");
}

void disp_mat(double **a, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++) {
	Rprintf("\n");
	for(j=0; j < m; j++) Rprintf("%10.8f ",a[i][j]);
    }
    Rprintf("\n");
}

void R_find_D_scale(double *rr, double *kkappa, double *ttau, int *llength,
		    double *sscale, double *cc, int *iipsi, int *ttype, double *rel_tol,
		    int *max_k, int *converged)
{
    /* compute D_scale using iterative algorithm
       type: 1: d1
       2: d2
       3: dt1
       4: dt2
    */
    *converged = 0;
    for (int k=0; k < *max_k; k++) {
	double scale = *sscale, tsum1 = 0, tsum2 = 0;
	// calculate weights
	for(int i=0; i < *llength; i++) {
	    double a, w = wgt(rr[i] / ttau[i] / scale, cc, *iipsi);
	    switch(*ttype) {
	    case 1: // d1
		a = rr[i]/ttau[i];
		tsum1 += a*a*w;
		tsum2 += w; break;
	    case 2: // d2
		a = rr[i]/ttau[i]*w;
		tsum1 += a*a;
		tsum2 += w*w; break;
	    default:
	    case 3: // dt1
		tsum1 += rr[i]*rr[i]*w;
		tsum2 += w*ttau[i]*ttau[i]; break;
	    case 4: // dt2
		a = rr[i]*w;
		tsum1 += a*a;
		a = ttau[i]*w;
		tsum2 += a*a; break;
	    };
	}

	*sscale = sqrt(tsum1 / tsum2 / *kkappa);

	// Rprintf("\n type = %d, scale = %10.8f \n", *ttype, *sscale);

	if (fabs(scale - *sscale) < *rel_tol * fmax2(*rel_tol, scale)) {
	    *converged = 1;
	    break;
	}
    }
}

/* specialized function calc_fitted */
/* calculates fitted values from simulation output array. */
/* this is used to process simulation output in the */
/* lmrob_simulation vignette */
void R_calc_fitted(double *XX, double *bbeta, double *RR, int *nn, int *pp, int *nnrep,
		   int *nnproc, int *nnerr)
{
    unsigned long A, B, C, D, E;
    A = (unsigned long)*nnerr; B = (unsigned long)*nnproc; C = (unsigned long)*nnrep;
    D = (unsigned long)*nn; E = (unsigned long)*pp;
    // calculate fitted values over errstr, procstr and replicates
    for(unsigned long a = 0; a < A; a++) { // errstr
	for(unsigned long b = 0; b < B; b++) { // procstr
	    for(unsigned long c = 0; c < C; c++) { // replicates
		// check for NAs
		if (!ISNA(bbeta[c + /* 0*C + */  b*C*E + a*B*E*C])) {
		    for(unsigned long d = 0; d < D; d++) { // observations
			RR[d + c*D + b*C*D + a*B*C*D] = 0; // initialize result
			for(unsigned long e = 0; e < E; e++) { // predictors
			    RR[d + c*D + b*C*D + a*B*C*D] += bbeta[c + e*C + b*C*E + a*B*E*C] *
				XX[d + e*D + c*E*D + a*C*E*D];
			}
		    }
		}
	    }
	}
    }
}
