/*
 *  Copyright (C) 2014	Martin Maechler, ETH Zurich
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#include <inttypes.h>
/*        ^^^^^^^^^^ is supposedly more common and standard than
 * #include <stdint.h>
 * or #include <sys/types.h> */
/* --> int64_t ; if people don't have the above, they can forget about it.. */
/* #include "int64.h" */

#include <Rmath.h> /* -> <math.h> and much more */

// Interface routines to be called via .C(), .Call() :
#include "robustbase.h"
//-> <R.h>, <Rinternals.h>  -> XLENGTH, R_xlen_t

/* Smooth Weighting Function -- typically for computing weights from large distances
 * -------------------------
 *                          \   <-- quartic polynomial here
 *                           \
 *                            ---------------------------
 * In fact a 2-parameter generalization of Tukey's 1-parameter "biweight"
 *
 * --- see also psi, rho, ... utilities  in ./lmrob.c
 */

double wgt_flex(double x, double c, double h) {
    double h2 = h/2.;
    x = fabs(x);
    if (x >= c+h2) return 0. ;
    if (x <= c-h2) return 1. ;
    // non-trivial {biweight like} down weighting:
    x = (x - (c-h2)) / h; // is in (0, 1)
    x = 1 - x*x;
    return x*x; // = (1 - ((|x| - (c - h/2))/ h)^2)^2  {in original 'x'}
}


SEXP R_wgt_flex(SEXP x_, SEXP c_, SEXP h_) { // TODO?: add , SEXP keep_attributes
    /*
     * Calculate Flexible weight function for vectorized x
     */
    int nprot = 1;
    if (isInteger(x_)) { x_ = PROTECT(coerceVector(x_, REALSXP)); nprot++; }
    if (isInteger(c_)) { c_ = PROTECT(coerceVector(c_, REALSXP)); nprot++; }
    if (isInteger(h_)) { h_ = PROTECT(coerceVector(h_, REALSXP)); nprot++; }
    if (!isReal(x_)) error(_("Argument '%s' must be numeric or integer"), "x");
    if (!isReal(c_) || LENGTH(c_) != 1) error(_("Argument '%s' must be numeric or integer of length 1"), "c");
    if (!isReal(h_) || LENGTH(h_) != 1) error(_("Argument '%s' must be numeric or integer of length 1"), "h");
    R_xlen_t i, n = XLENGTH(x_);
    SEXP res = PROTECT(allocVector(REALSXP, n)); // the result
    double *x = REAL(x_), *r = REAL(res), c = asReal(c_), h = asReal(h_);

    for(i = 0; i < n; i++)
	r[i] = ISNAN(x[i]) ? x[i] : wgt_flex(x[i], c, h);

    /* if(asLogical(keep_attributes)) { */
    // do the "no exception" version of copyMostAttrib() in ..R/src/main/attrib.c
    /* } */
    UNPROTECT(nprot);
    return res;
}

