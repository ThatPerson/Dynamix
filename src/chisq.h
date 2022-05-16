#ifndef CHISQ_H
#define CHISQ_H

#include "models/model.h"
#include "models/impact.h"

int opts_to_bcpars(Decimal *opts, struct BCParameters *pars, struct Model *m, struct Residue *resid, int *violations);

int bcpars_init(struct BCParameters *pars, Decimal slow, Decimal fast, struct Model *m);
int bcpars_update(struct BCParameters *pars, Decimal slow, Decimal fast, struct Model *m);
int bcpars_tdep(struct BCParameters *pars, struct BCParameters *npars, Decimal tfacts, Decimal tfactf);
Decimal back_calc(struct Residue *resid, struct Relaxation *relax, struct Model *m, int *violations, struct BCParameters *pars);
int bcpars_clean(struct BCParameters *pars);
Decimal optimize_chisq(Decimal *opts, struct Residue *resid, struct Model *m, unsigned int params);

int back_calculate(struct Residue *resid, struct Model *m, char *filename, unsigned int params);

void check_S2_violations(struct BCParameters *pars, int *violations);

#endif
