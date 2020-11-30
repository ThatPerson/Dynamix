#ifndef CHISQ_H
#define CHISQ_H

double back_calc(long double * opts, struct Residue * resid, struct Relaxation * relax, struct Model * m, int *violations);
double optimize_chisq(long double * opts, struct Residue * resid, struct Model * m, unsigned int params);
int back_calculate(long double * opts, struct Residue * resid, struct Model * m, char *filename, unsigned int params);

#include "chisq.c"
#endif
