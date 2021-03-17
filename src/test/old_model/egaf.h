#ifndef EGAF_H
#define EGAF_H

#ifndef GAF_H
#include "models/gaf.h"
#endif

#undef MOD_EGAF
#undef MOD_EGAFT
#define MOD_EGAF 			9
#define MOD_EGAFT 			10

Decimal EGAF_15NR1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal S2f);
Decimal EGAF_15NR2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal S2f);
Decimal EGAF_13CR1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal S2f);
Decimal EGAF_13CR2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal S2f);

#endif
