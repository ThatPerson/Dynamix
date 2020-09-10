#ifndef EGAF_H
#define EGAF_H

#ifndef GAF_H
#include "models/gaf.h"
#endif

#undef MOD_EGAF
#undef MOD_EGAFT
#define MOD_EGAF 			9
#define MOD_EGAFT 			10

double EGAF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double S2f);
double EGAF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double S2f);
double EGAF_13CR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double S2f);
double EGAF_13CR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double S2f);

#include "egaf.c"
#endif
