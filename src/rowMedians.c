/***************************************************************************

 Authors: Adopted from rowQuantiles.c by R. Gentleman.

 Copyright Henrik Bengtsson, 2007;  Martin Maechler, 2014;  History --> EOF
 **************************************************************************/
#include <Rdefines.h>

//  Public methods:

SEXP rowMedians_Real   (SEXP x, int nrow, int ncol, int narm, int hasna, int byrow);
SEXP rowMedians_Integer(SEXP x, int nrow, int ncol, int narm, int hasna, int byrow);

void C_rowMedians_Real   (double* x, double* res,
			  int nrow, int ncol, int narm, int hasna, int byrow);
void C_rowMedians_Integer(int*    x, double* res,
			  int nrow, int ncol, int narm, int hasna, int byrow);
/*
TEMPLATE rowMedians_<Integer|Real>(...):
- SEXP rowMedians_Real(...);
- SEXP rowMedians_Integer(...);
 */
#define METHOD rowMedians

#define X_TYPE 'i'
#include "rowMedians_TYPE-template.h"

#define X_TYPE 'r'
#include "rowMedians_TYPE-template.h"

#undef METHOD


/* TODO: implement: hasNA in {NA,TRUE,FALSE}; and = NA <==> code should *check*

   R code {for error message}: ../R/comedian.R  */
SEXP R_rowMedians(SEXP x, SEXP naRm, SEXP hasNA, SEXP byRow, SEXP keepNms) {

  // Argument checking and "C type coercion":
  if (!isMatrix(x))
    error("Argument 'x' must be a matrix.");

  int narm = asLogical(naRm); // error if it ain't
  if (narm != TRUE && narm != FALSE)
    error("Argument 'naRm' must be either TRUE or FALSE.");

  int hasna = asLogical(hasNA); // error if it ain't
  if (hasna == NA_INTEGER)
      hasna = TRUE;// <- for now; TODO ? become smarter and check

  int byrow = INTEGER(byRow)[0];
  int keepnms = asLogical(keepNms);

  /* Get dimensions of 'x'. */
  SEXP ans = PROTECT(getAttrib(x, R_DimSymbol));
  int nrow, ncol;
  if (byrow) { // rowMedians
    nrow = INTEGER(ans)[0];
    ncol = INTEGER(ans)[1];
  } else { // colMedians
    nrow = INTEGER(ans)[1];
    ncol = INTEGER(ans)[0];
  }

  if (isReal(x)) {
    ans = rowMedians_Real(x, nrow, ncol, narm, hasna, byrow);
  } else if (isInteger(x)) {
    ans = rowMedians_Integer(x, nrow, ncol, narm, hasna, byrow);
  } else {
    UNPROTECT(1);
    error("Argument 'x' must be numeric (integer or double).");
  }
  if(keepnms) {
      SEXP xDnms = getAttrib(x, R_DimNamesSymbol);
      if(xDnms != R_NilValue) {
	  PROTECT(xDnms);
	  setAttrib(ans, R_NamesSymbol,
		    duplicate(VECTOR_ELT(xDnms, byrow ? 0 : 1)));
	  UNPROTECT(1);
      }
  }
  UNPROTECT(1);
  return(ans);
} /* R_rowMedians() */


/***************************************************************************
 HISTORY:
 2014-12-09 [M.Maechler]
 o Copied to 'robustbase' CRAN package - to replace many apply(*., 2, median)
   NB: 'Biobase' also contains rowQ = general row/col Quantiles
 o argument checking all in C
 o add 'keepNms' argument {and do keep names by default!}

 2013-01-13 [HB]
 o Added argument 'byRow' to rowMedians() and dropped colMedians().
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
