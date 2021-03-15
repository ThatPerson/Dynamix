#ifndef GAF_H
#define GAF_H

#ifndef EMF_H
#include "models/emf.h"
#endif

#undef MOD_GAF
#undef MOD_GAFT
#define MOD_GAF 			7
#define MOD_GAFT 			8

Decimal GAF_Dipolar_R1(Decimal omega_obs, Decimal omega_neigh, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal D);
Decimal GAF_CSA_R2(Decimal omega, \
						 Decimal w1, \
						 Decimal wr, \
						 Decimal taus, \
						 Decimal S2s, \
						 Decimal tauf, \
						 Decimal S2f, \
						 Decimal D2, \
						 Decimal  (*J_SD)(\
						 	Decimal,\
						 	Decimal, \
						 	Decimal, \
						 	Decimal, \
						 	Decimal)\
						);
Decimal GAF_Dipolar_R2(Decimal omega_obs, Decimal omega_neigh, Decimal w1, Decimal wr, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal D);
int GAF_S2(Decimal sig[3], struct Orient ** A, struct Orient ** B, Decimal * S2[], int length,  int mode);
Decimal GAF_15NR1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf);
Decimal GAF_15NR2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf);
Decimal GAF_13CR1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf);
Decimal GAF_13CR2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf);

#endif
