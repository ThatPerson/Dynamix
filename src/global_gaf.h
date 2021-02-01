#ifndef GLOBAL_H
#define GLOBAL_H

double optimize_global_chisq(long double *opts, struct Residue *r, struct Model *m, unsigned int params);
void *run_global_iteration(void *input);
int calc_global_errors(struct Model *m);
int run_global(struct Model *m);

#include "global_gaf.c"
#endif
