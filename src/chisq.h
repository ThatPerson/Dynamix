#ifndef CHISQ_H
#define CHISQ_H

double back_calc(long double * opts, struct Residue * resid, struct Relaxation * relax, unsigned int model, int *violations);
double optimize_chisq(long double * opts, struct Residue * resid, unsigned int model, unsigned int or_variations, unsigned int params);
int back_calculate(long double * opts, struct Residue * resid, unsigned int model, unsigned int or_variations, char *filename, unsigned int params);

#include "chisq.c"
#endif
