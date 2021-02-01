#ifndef GLOBAL_H
#define GLOBAL_H

double optimize_global_chisq(long double *opts, struct Residue *r, struct Model *m, unsigned int params);
void *run_global(struct Model *m);

#include "global_gaf.c"
#endif
