#ifndef GLOBAL_H
#define GLOBAL_H

double optimize_global_chisq(double *opts, struct Residue *r, struct Model *m, unsigned int params);
int run_global_iteration(struct Model *m, int i);
int calc_global_errors(struct Model *m);
int run_global(struct Model *m);

#include "global_gaf.c"
#endif
