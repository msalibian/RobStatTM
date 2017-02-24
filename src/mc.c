/*
    Algorithm for the skewness estimator medcouple (MC)
    --------------------------------------------------
    ( originally matlabmc.c and also  mc/mcrsoft/spmc.c )
*/

#include <stdlib.h>
#include <math.h>

#include <inttypes.h>
// -> int64_t

#include <Rmath.h>
/* -> fmax2(.,.) */
#include <R_ext/Utils.h>


/* Interface routines to be called via .C() and those from API : */
#include "robustbase.h"
/*
  including

 whimed_i(a,iw,n): the weighted high median of an array a of length n,
		   using the positive integer weights iw[].

 * which is in ./wgt_himed.c_templ
 *		 ~~~~~~~~~~~~~~~~~
*/

/* Includes the auxiliary function

   h_kern(a,b, ai,bi,ab, eps):	the values h(a,b)  needed to compute the mc
*/
static
double h_kern(double a, double b, int ai, int bi, int ab, double eps);


void mc_C(double *z, int *in, double *eps, int *iter, double *out)
{
    *out = mc_C_d(z, *in, eps, iter);
    return;
}

/* MM:	The tolerance  'eps1' and 'eps2' can now be passed from R;
 *	the original code had only one 'eps' for both and hardcoded
 *	   eps =  0.0000000000001;  (== 1e-13 )
 *
 * MK:  eps1: for (relative) "equality" checks
 *      eps2: used to check for over- and underflow, respectively
 *      therefore I suggest eps1 = DBL_EPS and eps2 = DBL_MIN
 */
double mc_C_d(double *z, int n, double *eps, int *iter)
{
/* NOTE:
    eps	  = c(eps1, eps2)
    iter := c(maxit, trace.lev)  as input
          = c(it, converged)     as output
*/
    int trace_lev = iter[1], it = 0;
    Rboolean converged = TRUE;
    double medc; // "the" result
    static const double Large = DBL_MAX / 4.;

    if (n < 3) {
	medc = 0.; goto Finish;
    }
    /* copy data before sort()ing in place, also reflecting it -- dealing with +-Inf.
       NOTE: x[0] "empty" so we can use 1-indexing below */
    double *x  = (double *) R_alloc(n+1, sizeof(double));
    x[0] = 0;
    for (int i = 0; i < n; i++) {
	double zi = z[i];
	x[i+1] = - ((zi == R_PosInf) ? Large :
		    (zi == R_NegInf ? -Large : zi));
    }

    R_rsort(&x[1], n); /* full sort */

    double xmed; // := median( x[1:n] ) = - median( z[0:(n-1)] ):
    if (n%2) { /* n even */
	xmed = x[(n/2)+1];
    }
    else { /* n  odd */
	int ind = (n/2);
	xmed = (x[ind] + x[ind+1])/2;
    }

    if (fabs(x[1] - xmed) < eps[0] * (eps[0] + fabs(xmed))) {
	medc = -1.; goto Finish;
    } else if (fabs(x[n] - xmed) < eps[0] * (eps[0] + fabs(xmed))) {
	medc =	1.; goto Finish;
    }
    /* else : median is not at the border ------------------- */

    if(trace_lev)
	Rprintf("mc_C_d(z[1:%d], trace_lev=%d): Median = %g (not at the border)\n",
		n, trace_lev, -xmed);

    int i,j;
    /* center x[] wrt median --> such that then  median( x[1:n] ) == 0 */
    for (i = 1; i <= n; i++)
	x[i] -= xmed;

    /* Now scale to inside [-0.5, 0.5] and flip sign such that afterwards
    *  x[1] >= x[2] >= ... >= x[n] */
    double xden = -2 * fmax2(-x[1], x[n]);
    for (i = 1; i <= n; i++)
	x[i] /= xden;
    xmed /= xden;
    if(trace_lev >= 2)
	Rprintf(" x[] has been rescaled (* 1/s) with s = %g\n", -xden);

    j = 1;
    double x_eps = eps[0] * (eps[0] + fabs(xmed));
    while (j <= n && x[j] > x_eps) { /* test relative to xmed */
	/* x1[j] = x[j]; */
	j++;
    }
    if(trace_lev >= 2)
	Rprintf("   x1[] := {x | x_j > x_eps = %g}    has %d (='j-1') entries\n",
		x_eps, j-1);
    i = 1;
    double *x2 = x+j-1; /* pointer -- corresponding to  x2[i] = x[j]; */
    while (j <= n && x[j] > -x_eps) { /* test relative to xmed */
	/* x1[j] = x[j]; */
        /* x2[i] = x[j]; */
        j++;
        i++;
    }
    /* now  x1[] := {x | x_j > -eps}  also includes the median (0) */
    if(trace_lev >= 2)
        Rprintf("'median-x' {x | -eps < x_i <= eps} has %d (= 'k') entries\n",
		i-1);
    int h1 = j-1, /* == size of x1[] == the sum of those two sizes above */
    /* conceptually,  x2[] := {x | x_j <= eps}   (which includes the median 0) */
	h2 = i + (n-j);// == size of x2[] == maximal size of whimed() arrays

    if(trace_lev)
	Rprintf("  now allocating 2+5 work arrays of size (1+) h2=%d each:\n", h2);
    /* work arrays for whimed_i() :  allocate *once* only !! */
    double *acand  = (double *) R_alloc(h2, sizeof(double)),
	   *a_srt  = (double *) R_alloc(h2, sizeof(double));

    int    *iw_cand= (int *)	R_alloc(h2, sizeof(int)),
    /* work arrays for the fast-median-of-table algorithm:
     *  currently still with  1-indexing */
	*left  = (int *) R_alloc((h2+1), sizeof(int)),
	*right = (int *) R_alloc((h2+1), sizeof(int)),
	*p     = (int *) R_alloc((h2+1), sizeof(int)),
	*q     = (int *) R_alloc((h2+1), sizeof(int));

    for (i = 1; i <= h2; i++) {
	left [i] = 1;
	right[i] = h1;
    }
    int64_t nr = ((int64_t) h1) * ((int64_t) h2), // <-- careful to *NOT* overflow
	knew = nr/2 +1;
    if(trace_lev >= 2)
	Rprintf(" (h1,h2, nr, knew) = (%d,%d, %.0f, %.0f)\n",
		h1,h2, (double)nr, (double)knew);

    double trial = -2./* -Wall */;
    double *work= (double *) R_alloc(n, sizeof(double));
    int	   *iwt = (int *)    R_alloc(n, sizeof(int));
    Rboolean IsFound = FALSE;
    int nl = 0,
	neq = 0;
    /* MK:  'neq' counts the number of observations in the
     *      inside the tolerance range, i.e., where left > right + 1,
     *      since we would miss those when just using 'nl-nr'.
     *      This is to prevent index overflow in work[] later on.
     *      left might be larger than right + 1 since we are only
     *      testing with accuracy eps_trial and therefore there might
     *      be more than one observation in the `tolerance range`
     *      between < and <=.
     */
    while (!IsFound && (nr-nl+neq > n) && it < iter[0])
    {
	int64_t sum_p, sum_q;
	it++;
	j = 0;
	for (i = 1; i <= h2; i++)
	    if (left[i] <= right[i]) {
		iwt[j] = right[i] - left[i]+1;
		int k = left[i] + (iwt[j]/2);
		work[j] = h_kern(x[k], x2[i], k, i, h1+1, eps[1]);
		j++;
	    }
	if(trace_lev >= 4) {
	    Rprintf(" before whimed(): work and iwt, each [0:(%d-1)]:\n", j);
	    if(j >= 100) {
		for(i=0; i < 90; i++) Rprintf(" %8g", work[i]); Rprintf("\n  ... ");
		for(i=j-4; i < j; i++)Rprintf(" %8g", work[i]); Rprintf("\n");
		for(i=0; i < 90; i++) Rprintf(" %8d", iwt [i]); Rprintf("\n  ... ");
		for(i=j-4; i < j; i++)Rprintf(" %8d", iwt [i]); Rprintf("\n");
	    } else { // j <= 99
		for(i=0; i < j; i++) Rprintf(" %8g", work[i]); Rprintf("\n");
		for(i=0; i < j; i++) Rprintf(" %8d", iwt [i]); Rprintf("\n");
	    }
	}
	trial = whimed_i(work, iwt, j, acand, a_srt, iw_cand);
	double eps_trial = eps[0] * (eps[0] + fabs(trial));
	if(trace_lev >= 3)
	    Rprintf("%2s it=%2d, whimed(*, n=%6d)=%11.5g ", " ", it, j, trial);

	j = 1;
	for (i = h2; i >= 1; i--) {
	    while (j <= h1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1]) - trial > eps_trial) {
		// while (j <= h1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1]) > trial) {
		if (trace_lev >= 5)
		    Rprintf("\nj=%3d, i=%3d, x[j]=%g, x2[i]=%g, h=%g",
			    j, i, x[j], x2[i],
			    h_kern(x[j],x2[i],j,i,h1+1,eps[1]));
		j++;
	    }
/* 	    for(; j <= h1; j++) { */
/* 		register double h = h_kern(x[j],x2[i],j,i,h1+1,eps[1]); */
/* 		if(h > trial) break; */
/* 	    } */
	    p[i] = j-1;
	}
	j = h1;
	for (i = 1, sum_p=0, sum_q=0; i <= h2; i++) {
	    while (j >= 1 && trial - h_kern(x[j],x2[i],j,i,h1+1,eps[1]) > eps_trial)
		// while (j >= 1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1]) < trial)
		j--;
	    q[i] = j+1;

	    sum_p += p[i];
	    sum_q += j;/* = q[i]-1 */
	}

	if(trace_lev >= 3) {
	    if (trace_lev == 3)
		Rprintf("sum_(p,q)= (%.0f,%.0f)", (double)sum_p, (double)sum_q);
	    else { /* trace_lev >= 4 */
		Rprintf("\n%3s p[1:%d]:", "", h2);
		Rboolean lrg = h2 >= 100;
		int i_m = lrg ? 95 : h2;
		for(i = 1; i <= i_m; i++) Rprintf(" %2d", p[i]); if(lrg) Rprintf(" ...");
		Rprintf(" sum=%4.0f\n%3s q[1:%d]:", (double)sum_p, "", h2);
		for(i = 1; i <= i_m; i++) Rprintf(" %2d", q[i]); if(lrg) Rprintf(" ...");
		Rprintf(" sum=%4.0f\n", (double)sum_q);
	    }
	}

	if (knew <= sum_p) {
	    if(trace_lev >= 3)
		Rprintf("; sum_p >= kn\n");
	    for (i = 1, neq = 0; i <= h2; i++) {
		right[i] = p[i];
		if (left[i] > right[i]+1) neq += left[i]-right[i]-1;
	    }
	    nr = sum_p;
	}
	else { /* knew > sum_p */
	    IsFound = (knew <= sum_q); /* i.e. sum_p < knew <= sum_q */;

	    if(trace_lev >= 3)
		Rprintf("; s_p < kn ?<=? s_q: %s\n", IsFound ? "TRUE": "no");
	    if(IsFound) {
		medc = trial;
	    } else { /*	 knew > sum_q */
	        for (i = 1; i <= h2; i++) {
		    left[i] = q[i];
		    if (left[i] > right[i]+1) neq += left[i]-right[i]-1;
		}
		nl = sum_q;
	    }
	}
	R_CheckUserInterrupt();

    } /* end while loop */

    converged = IsFound || (nr-nl+neq <= n);
    if(!converged) {
	warning("maximal number of iterations (%d =? %d) reached prematurely\n",
		 iter[0], it);
	/* still: */
	medc = trial;
    }

    if (converged && !IsFound) { /* e.g., for  mc(1:4) : */
	j = 0;
	for (i = 1; i <= h2; i++) {
	    if (left[i] <= right[i]) {
		for (int k = left[i]; k <= right[i]; k++) {
		    work[j] = -h_kern(x[k],x2[i],k,i,h1+1,eps[1]);
		    j++;
		}
	    }
	}
	if(trace_lev)
	    Rprintf("  not found [it=%d,  (nr,nl) = (%d,%d)],"
		    " -> (knew-nl, j) = (%d,%d)\n",
		    it, nr, nl, knew-nl, j);
	/* using rPsort(work, n,k), since we don't need work[] anymore:*/
	rPsort(work, /* n = */ j, /* k = */ knew-nl-1);
	medc = - work[knew-nl-1];
    }

    if(converged && trace_lev >= 2)
	Rprintf("converged in %d iterations\n", it);

Finish:
    iter[0] = it; /* to return */
    iter[1] = converged;

    return medc;

} /* end{ mc_C_d } */


static
double h_kern(double a, double b, int ai, int bi, int ab, double eps)
{
/*     if (fabs(a-b) <= DBL_MIN) */
    /* check for zero division and positive b */
    if (fabs(a-b) < 2.0*eps || b > 0)
	return sign((double)(ab - (ai+bi)));

    /* else */
    return (a+b)/(a-b);
}


/* Local variables section

 * Local variables:
 * mode: c
 * kept-old-versions: 12
 * kept-new-versions: 20
 * End:
 */
