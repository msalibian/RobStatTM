/* included 2 x from ./rowMedians.c
                       ~~~~~~~~~~~~
 ***********************************************************************
 TEMPLATE:
  SEXP rowMedians_<Integer|Real>(...)

 GENERATES:
  SEXP rowMedians_Integer(SEXP x, int nrow, int ncol, int narm, int hasna, int byrow)
  SEXP rowMedians_Real(SEXP x, int nrow, int ncol, int narm, int hasna, int byrow)

 Arguments:
   The following macros ("arguments") should be defined for the
   template to work as intended.

  - METHOD: the name of the resulting function
  - X_TYPE: 'i' or 'r'

 Authors:
  Adopted from rowQuantiles.c by R. Gentleman.
  Template by Henrik Bengtsson.

 Copyright: Henrik Bengtsson, 2007-2013;  Martin Maechler 2014
 ***********************************************************************/
#include <R.h>
// MM{FIXME}: only <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

/* Expand arguments:
    X_TYPE => (X_C_TYPE, X_IN_C, X_ISNAN, [METHOD_NAME])
 */
#include "templates-types.h"

#if X_TYPE == 'i'
  #define PSORT iPsort
#elif X_TYPE == 'r'
  #define PSORT rPsort
#endif

SEXP METHOD_NAME(SEXP x, int nrow, int ncol, int narm, int hasna, int byrow)
{
  /* R allocate a double vector of length 'nrow'
     Note that 'nrow' means 'ncol' if byrow=FALSE. */
  SEXP ans = PROTECT(allocVector(REALSXP, nrow));

  C_METHOD_NAME(X_IN_C(x), REAL(ans), nrow, ncol, narm, hasna, byrow);

  UNPROTECT(1);
  return(ans);
}

// "C-only" method [SEXP free, can be called from "pure C" :
void C_METHOD_NAME(X_C_TYPE *x, double *res,
		   int nrow, int ncol, int narm, int hasna, int byrow)
{
  Rboolean isOdd;
  int ii, jj, kk, qq;
  int *colOffset;
  X_C_TYPE value,
      /* R allocate memory for the 'rowData'.  This will be
	 taken care of by the R garbage collector later on. */
      *rowData = (X_C_TYPE *) R_alloc(ncol, sizeof(X_C_TYPE));

  if (!hasna) // If there are no missing values, don't try to remove them
    narm = FALSE;

  /* When narm is false, isOdd and qq are the same for all rows */
  if (!narm) {
    isOdd = (ncol % 2 == 1);
    qq = (int)(ncol/2) - 1;
  } else {
    isOdd = FALSE;
    qq = 0;
  }

  value = 0;

  /* Pre-calculate the column offsets */
  colOffset = (int *) R_alloc(ncol, sizeof(int));

  // HJ begin
  if (byrow) {
    for(jj=0; jj < ncol; jj++)
      colOffset[jj] = (int)jj*nrow;
  } else {
    for(jj=0; jj < ncol; jj++)
      colOffset[jj] = jj;
  }
  // HJ end

  if (hasna) {
    for(ii=0; ii < nrow; ii++) {
      if(ii % 1000 == 0)
        R_CheckUserInterrupt();

      int rowIdx = byrow ? ii : ncol*ii; //HJ

      kk = 0;  /* The index of the last non-NA value detected */
      for(jj=0; jj < ncol; jj++) {
        value = x[rowIdx+colOffset[jj]]; //HJ

        if (X_ISNAN(value)) {
          if (!narm) {
            kk = -1;
            break;
          }
        } else {
          rowData[kk] = value;
          kk = kk + 1;
        }
      }

      if (kk == 0) {
        res[ii] = R_NaN;
      } else if (kk == -1) {
        res[ii] = R_NaReal;
      } else {
        /* When narm is true, isOdd and qq may change with row */
        if (narm) {
          isOdd = (kk % 2 == 1);
          qq = (int)(kk/2) - 1;
        }

        /* Permute x[0:kk-1] so that x[qq] is in the correct
           place with smaller values to the left, ... */
        PSORT(rowData, kk, qq+1);
        value = rowData[qq+1];

        if (isOdd) {
          res[ii] = (double)value;
        } else {
          if (narm || !X_ISNAN(value)) {
            /* Permute x[0:qq-2] so that x[qq-1] is in the correct
               place with smaller values to the left, ... */
            PSORT(rowData, qq+1, qq);
            if (X_ISNAN(rowData[qq]))
              res[ii] = R_NaReal;
            else
              res[ii] = ((double)(rowData[qq] + value))/2;
          } else {
            res[ii] = (double)value;
          }
        }
      }
    } // for(..)
  } else { // no NAs
    for(ii=0; ii < nrow; ii++) {
      if(ii % 1000 == 0)
        R_CheckUserInterrupt();

      int rowIdx = byrow ? ii : ncol*ii; //HJ

      for(jj=0; jj < ncol; jj++)
        rowData[jj] = x[rowIdx+colOffset[jj]]; //HJ

      /* Permute x[0:ncol-1] so that x[qq] is in the correct
         place with smaller values to the left, ... */
      PSORT(rowData, ncol, qq+1);
      value = rowData[qq+1];

      if (isOdd) {
        res[ii] = (double)value;
      } else {
        /* Permute x[0:qq-2] so that x[qq-1] is in the correct
           place with smaller values to the left, ... */
        PSORT(rowData, qq+1, qq);
        res[ii] = (double)((rowData[qq] + value))/2;
      }
    } // for(..)
  } /* if (hasna ...) */
}

/* Undo template macros */
#undef PSORT
#include "templates-types_undef.h"


/***************************************************************************
 HISTORY:
 2014-12-09 [MMaechler]
 o do not use '== TRUE' '== FALSE' -- as we have no NA here
 o resolve REAL(ans) outside for(ii ..)
 o add "SEXP-free" C routines, others can call: C_rowMedians_(Real|Integer)

 2013-04-23 [HB]
  o BUG FIX: The integer template of rowMedians_<Integer|Real>() would
    not handle ties properly.  This was because ties were calculated as
    '(double)((rowData[qq] + value)/2)' instead of
    '((double)(rowData[qq] + value))/2'.
 2013-01-13 [HB]
  o Merged rowMedians_Integer() and rowMedians_Read() into template
    rowMedians_<Integer|Real>().
 2013-01-13 [HB]
 o Using internal arguments 'by_row' instead of 'by_column'.
 2011-12-11 [HB]
 o BUG FIX: rowMediansReal(..., na.rm=TRUE) did not handle NaN:s, only NA:s.
   Note that NaN:s does not exist for integers.
 2011-10-12 [HJ]
 o Added colMedians().
 o Now rowMediansInteger/Real() can operate also by columns, cf. argument
   'by_column'.
 2007-08-14 [HB]
 o Added checks for user interrupts every 1000 line.
 o Added argument 'hasNA' to rowMedians().
 2005-12-07 [HB]
 o BUG FIX: When calculating the median of an even number (non-NA) values,
    the length of the second sort was one element too short, which made the
    method to freeze, i.e. rPsort(rowData, qq, qq) is now (...qq+1, qq).
 2005-11-24 [HB]
  o By implementing a special version for integers, there is no need to
    coerce to double in R, which would take up twice the amount of memory.
  o rowMedians() now handles NAs too.
  o Adopted from rowQuantiles.c in Biobase of Bioconductor.
 **************************************************************************/
