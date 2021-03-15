#ifndef CHISQ_H
#define CHISQ_H

Decimal back_calc(Decimal * opts, struct Residue * resid, struct Relaxation * relax, struct Model * m, int *violations);
Decimal optimize_chisq(Decimal * opts, struct Residue * resid, struct Model * m, unsigned int params);
int back_calculate(Decimal * opts, struct Residue * resid, struct Model * m, char *filename, unsigned int params);

#endif
