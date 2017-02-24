#include <R.h>

void F77_SUB(rndstart)(void)  { GetRNGstate(); }
void F77_SUB(rndend)(void)    { PutRNGstate(); }
double F77_SUB(normrnd)(void) { return norm_rand(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }
