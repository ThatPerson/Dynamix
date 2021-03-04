#ifndef CHISQ_H
#define CHISQ_H

double back_calc(double * opts, struct Residue * resid, struct Relaxation * relax, struct Model * m, int *violations);
double optimize_chisq(double * opts, struct Residue * resid, struct Model * m, unsigned int params);
int back_calculate(double * opts, struct Residue * resid, struct Model * m, char *filename, unsigned int params);

#include "chisq.c"
#endif
