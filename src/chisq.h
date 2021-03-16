#ifndef CHISQ_H
#define CHISQ_H

#include "models/model.h"

int opts_to_bcpars(Decimal *opts, struct BCParameters *pars, unsigned int model, struct Residue *resid, int *violations);
int bcpars_init(struct BCParameters *pars, Decimal slow, Decimal fast);
Decimal back_calc(struct Residue * resid, struct Relaxation * relax, struct Model * m, int *violations, struct BCParameters *pars);
Decimal optimize_chisq(Decimal * opts, struct Residue * resid, struct Model * m, unsigned int params);
int back_calculate(Decimal * opts, struct Residue * resid, struct Model * m, char *filename, unsigned int params);

#endif
