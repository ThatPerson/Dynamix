#ifndef CHISQ_H
#define CHISQ_H

double optimize_chisq(long double * opts, struct Residue * resid, unsigned int model);
int back_calculate(long double * opts, struct Residue * resid, unsigned int model, char *filename);

#include "chisq.c"
#endif
