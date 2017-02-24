
#include <R_ext/Rdynload.h>
#include "robustbase.h"


#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static R_NativePrimitiveArgType Qn0_t[] = {
    REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType Sn0_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType mc_C_t[] = {
    REALSXP, INTSXP, REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType wgt_himed_i_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP
};
static R_NativePrimitiveArgType wgt_himed_t[] = {
    REALSXP, INTSXP, REALSXP, REALSXP
};


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

static const R_CMethodDef CEntries[]  = {
    CDEF(Qn0),
    CDEF(Sn0),
    CDEF(mc_C),
    CDEF(wgt_himed_i),
    CDEF(wgt_himed),
    CDEF(R_lmrob_S),
    CDEF(R_lmrob_MM),
    CDEF(R_find_D_scale),
    CDEF(R_calc_fitted),
    CDEF(R_lmrob_M_S),
    CDEF(R_subsample),
    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {
    CALLDEF(R_rho_inf, 2), // -> lmrob.c
    CALLDEF(R_psifun, 4),
    CALLDEF(R_chifun, 4),
    CALLDEF(R_wgtfun, 3),
    CALLDEF(R_wgt_flex, 3), // -> rob-utils.c
    CALLDEF(R_rowMedians, 5),// -> rowMedians.c [Biobase also has rowQ for quantiles]
    {NULL, NULL, 0}
};


static R_FortranMethodDef FortEntries[] = {
    {"rffastmcd", (DL_FUNC) &F77_SUB(rffastmcd), 49},/* ./rffastmcd.f */
    {"rfltsreg",  (DL_FUNC) &F77_SUB(rfltsreg), 41}, /* ./rfltsreg.f */
    {"rllarsbi",  (DL_FUNC) &F77_SUB(rllarsbi), 18}, /* ./rllarsbi.f */
    {NULL, NULL, 0}
};

void R_init_robustbase(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
