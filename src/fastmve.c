

/* FILE: fast-mve.c
//
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #include "R.h"
// #include "R_ext/Linpack.h"
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>


// #include "RobStatTM.h"

#define INFI 1e+20


void r_fast_mve(double *xx, int *nn, int *pp, int *nnsamp,
			int *nsingular, double *ctr, double *ccov, double *scale,
			int *best_ind, int *nnind, int *nn2)
{

void r_find_k_smallest(double *a, int n, int k, int *ind, double *tmp);
void r_mean_cov_mah_sample(double *x, int *n, int *p,
		int *indices, int *nind, double *xw, double *mean,
		double *cov, double *mah_d, double *det,
		int *pivot, double *qraux, double *work, int *rank,
		int *compute_cov, int *compute_distances, int *compute_det);
double median(double *x, int n, double *aux);
void sample_noreplace_mve(int *x, int n, int k, int *all_ind);

int *indi, iter, *pivot, rank, *all_ind;
int compute_cov, compute_det, compute_distances;
register int i, j;
int n = *nn, p = *pp, nsamp = *nnsamp, n2 = *nn2, nind = *nnind;
double s0 = INFI, s, *best_ctr, det1, det2;
double  *di, *tmp, *qraux, *tmp2; 
double *cov_vec;

all_ind  = (int *) malloc( n * sizeof(int) );
cov_vec  = (double*) calloc( p * p, sizeof(double) );
qraux    = (double *) malloc( p * sizeof(double) );
tmp2     = (double *) malloc( 2 * p * sizeof(double) ); /* work */
pivot    = (int *) malloc( p * sizeof(int) );
indi     = (int *) calloc( n , sizeof(int) );
best_ctr = (double *) malloc( p * sizeof(double) );
di       =  (double *) malloc( n * sizeof(double) );
tmp      =  (double *) malloc( n * p * sizeof(double) );

GetRNGstate(); /* set the seed from R? */

for(iter=0; iter < nsamp; iter++) {
	compute_distances = 1;
	compute_det = 0; 
	compute_cov = 0;
	R_CheckUserInterrupt();	
	rank = 0;
	sample_noreplace_mve(indi, n, nind, all_ind); 
	r_mean_cov_mah_sample(xx, nn, pp, indi, &nind, tmp, ctr, 
			cov_vec, di, &det1, pivot, qraux, tmp2, &rank,
			&compute_cov, &compute_distances, &compute_det);
	if( rank == p )  { 
		r_find_k_smallest(di, n, n2, indi, tmp);
		compute_distances = 1;
		compute_det = 1;
		compute_cov = 1;
		r_mean_cov_mah_sample(xx, nn, pp, indi, &n2, tmp, ctr, 
				cov_vec, di, &det1, pivot, qraux, tmp2, &rank,
				&compute_cov, &compute_distances, &compute_det);
		if( rank == p ) {
			det1 = det1 * det1 / pow( (double) n2 - 1.0, (double) p);
			/* det1 = |cov matrix| */
			det2 = pow(det1, 1.0 / (double) p );
			s = median(di, n, tmp) * det2;
			if( s < s0 )  {
				s0 = s;
				for(i=0; i<p; i++) {
					best_ctr[i] = ctr[i];
					for(j=0; j<p; j++)
						ccov[ j*p + i] = cov_vec[i + j*p] / det2;
				};
				for(i=0;i<n2;i++) best_ind[i] = indi[i]+1;
			};
		};
	} else *nsingular = *nsingular + 1;
}; /* end iter = nsamp */

for(i=0;i<p;i++) ctr[i] = best_ctr[i];
*scale = s0;

free(all_ind); free(qraux); free(pivot); free(tmp2); 
free(cov_vec); free(indi); free(di); free(tmp); free(best_ctr);
}

// double median(double *x, int n, double *aux) 
// {
// double kthplace(double *,int,int);
// double t;
// register int i;
// for(i=0;i<n;i++) aux[i]=x[i];
// if ( (n/2) == (double) n / 2 )
// 	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2 ;
// else	t = kthplace(aux,n, n/2+1 ) ;
// return(t);
// }

/* double kthplace(double *a, int n, int k)
{
int jnc,j;
int l,lr;
double ax,w;
k--;
l=0;
lr=n-1;
while (l<lr)
	{ ax=a[k];
	  jnc=l;
	  j=lr;
	  while (jnc<=j)
		{ while (a[jnc] < ax) jnc++;
		  while (a[j] > ax) j--;
		  if (jnc <= j)
			{ w=a[jnc];
			  a[jnc]=a[j];
			  a[j]=w;
			  jnc++;
			  j--;
			};
		};
	  if (j<k) l=jnc;
	if (k<jnc) lr=j;
	};
return(a[k]);
}
*/

void r_mean_cov_mah_sample(double *x, 
		int *n, int *p,
		int *indices, int *nind, double *xw, double *mean,
		double *cov, double *mah_d, double *det,
		int *pivot, double *qraux, double *work, int *rank,
		int *compute_cov, int *compute_distances,
		int *compute_det)
{
/* need double(p) in qraux, double (2*p) in work,
 * int(p) in pivot, double(p*p) in cov, double(n*p) in xw
 * int( *nind) in indices, double(p) in mean, double(n) in mah_d
 * det = sqrt( det( t(x[indices,]) %*% x[indices,] ) )
 */
double r_mah(double *xr, int nnew, int p, double *x, double *work);
int i, j, k;
int nn = *n, pp = *p, nnind = *nind;
double tol = 1e-7, s=0;

/* compute the mean, put submatrix into xw, center its columns */
for(j=0;j<pp;j++) {
	mean[j] = 0.0;
	for(i=0;i<nnind;i++) 
		mean[j] += (xw[i + j*nnind] = x[ indices[i] + j*nn ]) / (double) nnind;
	for(i=0;i<nnind;i++) 
		xw[i + j*nnind] -= mean[j];
};

/* QR decomposition of the submatrix */
F77_CALL(dqrdc2)(xw, nind, nind, p, &tol, rank, qraux, pivot, work);

/* build the cov matrix of the subsample using the QR decomp */
if(*compute_cov) {
for(i=0;i<pp;i++)
	for(j=i;j<pp;j++) {
		s = 0;
		for(k=0;k<=i;k++) s += xw[k + j*nnind] * xw[k + i*nnind];
		cov[j + i*pp] = cov[i + j*pp] = (s / (double) (nnind - 1));
	};
};
/* if full rank, compute det of cov matrix and mah distances of
 * all points based on mean and cov matrix of subsample */
/* det^2 = (nnind-1)^p * det( cov matrix ) */
if( *rank == pp) { 
	if(*compute_det) {
		*det = 1.0;
		for(j=0;j<pp;j++)
			*det *= fabs(xw[j + nnind*j]);
	};
	if(*compute_distances) {
		for(i=0;i<nn;i++) {
			for(j=0;j<pp;j++)
				qraux[j] = x[i + nn*j] - mean[j];
			mah_d[i] = r_mah(xw, nnind, pp, qraux, work);
		};
	};
} else
	*det = 0.0; 
}

void r_find_k_smallest(double *a, int n, int k, int *ind, double *tmp)
{
double kthplace(double *,int,int);
double aux;
int i,j;
for(i=0; i<n; i++) tmp[i] = a[i]; 
aux = kthplace(tmp, n, k);
j = 0;
for(i=0; i<n; i++) 
	if(a[i] <= aux) ind[j++] = i;
}

/* find the squared Mahalanobis distance to x via QR decomposition in xr. */
/* (c) B. Ripley */
double r_mah(double *xr, int nnew, int p, double *x, double *work)
{
    int i, j;
    double s, ss = 0.0;
    for(j = 0; j < p; j++) {
	s = x[j];
	if(j > 0) for(i = 0; i < j; i++) s -= work[i] * xr[i + nnew*j];
	work[j] = s / xr[j + nnew*j];
	ss += work[j] * work[j];
    }
    return(ss*(nnew-1));
}

/* Sampling k from 0:n-1 without replacement.
   (c) B. Ripley */
void sample_noreplace_mve(int *x, int n, int k, int *ind)
{
    int i, j, nn=n;

    for (i = 0; i < n; i++) ind[i] = i; 
    for (i = 0; i < k; i++) {
	j = nn * unif_rand(); 
	x[i] = ind[j];
	ind[j] = ind[--nn];
    }
}


