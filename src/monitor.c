// Flexible printing of informative messages from Fortran
// -----------------                         ------------
// other MM versions: ~/R/Pkgs/cobs99/src/monitor.c ~/R/Pkgs/lokern/src/monitor.c

#include <R.h>

/* called for  trace >= 2 : ----------------------------------------------- */

void F77_SUB(pr1mcd)(int *i_trace, int *n, int *nvar, int *nhallf, int *krep,
		     int *nmini, int *kmini)
{
    Rprintf("rffastmcd(n=%d, nvar=%d, nhallf=%d, krep=%d, nmini=%d, kmini=%d, i_trace=%d)\n",
	    *n, *nvar, *nhallf, *krep, *nmini, *kmini,
	    *i_trace);
}

void F77_SUB(pr2mcd)(Rboolean *part, Rboolean *all, // <- logical
		     int *kstep, int *ngroup, int *minigr, int *nhalf, int *nrep)
{
    Rprintf("pr[2]: (part=%d, all=%d); (kstep=%d, ngroup=%d, minigr=%d, nhalf=%d, nrep=%d)\n",
	    *part, *all, *kstep, *ngroup, *minigr, *nhalf, *nrep);
}


void F77_SUB(pr3mcd)(Rboolean *part, Rboolean *fine, int *final, // <- logical
		     int *nrep, int *nn, int *nsel, int *nhalf, int *kstep,
		     int *nmini, int *nmaxi)
{
    char* phase_kind = (*part)
	? ((*fine && !*final)
	   ? "fine (2 of 3)"
	   : ((*final) ? "final (3 of 3)" : "first (of 3)"))
	: ((*final) ? "final" : "one");
    Rprintf(" Main loop, phase[%s]:\n (nrep=%4d, nn=%4d, nsel=%4d, nhalf=%4d, kstep=%d, nmini=%d, nmaxi=%d)\n",
	    phase_kind, *nrep, *nn, *nsel, *nhalf,
	    *kstep, *nmini, *nmaxi);
}

void F77_SUB(prp1mcd)(int *n, int *ngroup, int *minigr, int *nhalf, int *nrep,
		      int mini[])
{
    // int mini[*kmini];
    Rprintf(" Partitioning n=%d into at most kmini groups: ngroup=%d, minigr=%d, nhalf=%d, nrep=%d;"
	    "\n groups are of sizes (",
	    *n, *ngroup, *minigr, *nhalf, *nrep);
    for(int j=0; j < *ngroup; j++) Rprintf(" %d", mini[j]);
    Rprintf(")\n");
}

void F77_SUB(pr9mcd)(int *ntot) {
    Rprintf(" -- finishing: total times = %d\n", *ntot);
}


/* called for  trace >= 3 : ----------------------------------------------- */

void F77_SUB(prgrmcd)(int *ii, int *nn, int *i_trace)
{
    Rprintf(" group ii = %d (nn = %d)%s\n", *ii, *nn,
	    (*i_trace >= 4) ? ": i=1..nrep loop: " : "");
}

void F77_SUB(pr4mcd)(int *i) { Rprintf(" i = %d "); }

void F77_SUB(pr5mcd)(int *step, int *ntot) {
    Rprintf("(step %d, tot=%d)", *step, *ntot);
}
