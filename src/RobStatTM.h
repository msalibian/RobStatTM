/* Copied from robustbase.h */

#include <R.h>
#include <Rinternals.h>
#include <complex.h>

/* For internationalized messages */
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP : String)
#endif

double complex erfz(double complex z);
double complex erfi(double complex z);
SEXP R_erfi(SEXP x);

SEXP R_rho_inf(SEXP cc, SEXP ipsi);

void R_lmrob_S(double *X, double *y, int *n, int *P,
	       int *nRes, double *scale, double *beta_s,
	       double *C, int *iipsi, double *bb,
	       int *best_r, int *Groups, int *N_group,
	       int *K_s, int *max_k, int *max_it_scale,
	       double *rel_tol, double *inv_tol,
               //     ^^^^^^^^^ = refine.tol in R
	       int* converged, int *trace_lev, int *mts, int *ss, int *cutoff);

void R_lmrob_M_S(double *X1, double *X2, double *y, double *res,
		 int *n, int *p1, int *p2, int *nRes, int *max_it_scale,
		 double *scale, double *b1, double *b2,
		 double *rho_c, int *ipsi, double *bb,
		 int *K_m_s, int *max_k, double *rel_tol, double *inv_tol,
		 int *converged, int *trace_lev,
		 int *orthogonalize, int *subsample,
		 int *descent, int *mts, int *ss);

void R_lmrob_MM(double *X, double *y, int *n, int *P,
		double *beta_initial, double *scale,
		double *beta_m, double *resid,
		int *max_it,
		double *rho_c, int *ipsi, double *loss, double *rel_tol,
		int *converged, int *trace_lev, int *mts, int *ss);

void R_subsample(const double *x, const double *y, int *n, int *m,
		 double *beta, int *ind_space, int *idc, int *idr,
		 double *lu, double *v, int *p,
		 double *_Dr, double *_Dc, int *_rowequ, int *_colequ,
		 int *status, int *sample, int *mts, int *ss, double *tol_inv,
		 int *solve);

SEXP R_psifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_);
SEXP R_chifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_);
SEXP R_wgtfun(SEXP x_, SEXP c_, SEXP ipsi_);


void R_find_D_scale(double *rr, double *kkappa, double *ttau, int *llength,
		    double *sscale, double *cc, int *iipsi, int *ttype, double *rel_tol,
		    int *max_k, int *converged);

void R_calc_fitted(double *XX, double *bbeta, double *RR, int *nn, int *pp, int *nnrep,
		   int *nnproc, int *nnerr);

void  F77_NAME(rllarsbi)(
    double *X, double *Y, int *N, int *NP, int *MDX, int *MDT,
    double *TOL, int *NIT, int *K, int *KODE, double *SIGMA, double *THETA,
    double *RS, double *SC1, double *SC2, double *SC3, double *SC4,
    double *BET0);

void r_fast_mve(double *xx, int *nn, int *pp, int *nnsamp,
                int *nsingular, double *ctr, double *ccov, double *scale,
                int *best_ind, int *nnind, int *nn2);

void F77_NAME(dqrdc2)(double *xw, int *nind, int *nind2, int *p,
              double *tol, int *rank, double *qraux, int *pivot, double *work);

