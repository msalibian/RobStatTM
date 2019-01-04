#include <R_ext/Rdynload.h>
#include "RobStatTM.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static R_NativePrimitiveArgType R_lmrob_S_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    /* rrhoc */ REALSXP, INTSXP, REALSXP,
    /* best_r */ INTSXP, INTSXP, INTSXP,
    /* K_s */ INTSXP, INTSXP, INTSXP,
    /* rel_tol*/ REALSXP, REALSXP,
    /* converged */ LGLSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_lmrob_MM_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP,
    /* beta_initial */ REALSXP, REALSXP,
    /* beta_m */ REALSXP, REALSXP,
    /* max_it */ INTSXP, REALSXP, INTSXP,
    /* loss */ REALSXP, REALSXP, LGLSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_find_D_scale_t[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    /* c */ REALSXP, INTSXP, INTSXP, REALSXP,
    /* max_k */ INTSXP, LGLSXP
};

static R_NativePrimitiveArgType R_calc_fitted_t[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_lmrob_M_S_t[] = {
    REALSXP, REALSXP, REALSXP, REALSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, REALSXP,
    INTSXP, INTSXP, REALSXP, REALSXP,
    LGLSXP, INTSXP,
    LGLSXP, LGLSXP, LGLSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_subsample_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP,
    REALSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, INTSXP, INTSXP,
    INTSXP, LGLSXP, INTSXP, INTSXP, REALSXP,
    LGLSXP
};


static R_NativePrimitiveArgType r_fast_mve_t[] = {
  REALSXP, INTSXP, INTSXP, INTSXP, 
  INTSXP, REALSXP, REALSXP, REALSXP, 
  INTSXP, INTSXP, INTSXP
};

static const R_CMethodDef CEntries[]  = {
    CDEF(r_fast_mve),
    CDEF(R_lmrob_S),
    CDEF(R_lmrob_MM),
    CDEF(R_find_D_scale),
    CDEF(R_calc_fitted),
    CDEF(R_lmrob_M_S),
    CDEF(R_subsample),
    {NULL, NULL, 0}
};


static const R_CallMethodDef CallEntries[] = {
    CALLDEF(R_rho_inf, 2),
    CALLDEF(R_psifun, 4),
    CALLDEF(R_chifun, 4),
    CALLDEF(R_wgtfun, 3),
    CALLDEF(R_erfi, 1),
    {NULL, NULL, 0}
};


static const R_FortranMethodDef FortEntries[] = {
    {"rslarsbi",  (DL_FUNC) &F77_SUB(rllarsbi), 18},
    {"dqrdc2", (DL_FUNC) &F77_SUB(dqrdc2), 9},
    {NULL, NULL, 0}
};


void R_init_RobStatTM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE); 
}
