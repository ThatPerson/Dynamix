#ifndef AIMF_H
#define AIMF_H

//#ifndef EMF_H
//#include "models/emf.h"
//#endif

#undef MOD_AIMF
#undef MOD_AIMFT
#define MOD_AIMF 			23
#define MOD_AIMFT 			24

Decimal AIMF_Dipolar_R1(Decimal omega_obs, Decimal omega_neigh, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal D);
Decimal AIMF_CSA_R2(Decimal omega, \
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
Decimal AIMF_Dipolar_R2(Decimal omega_obs, Decimal omega_neigh, Decimal w1, Decimal wr, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal D);
Decimal AIMF_15NR1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf);
Decimal AIMF_15NR2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf);
Decimal AIMF_13CR1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf);
Decimal AIMF_13CR2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf);

#endif
