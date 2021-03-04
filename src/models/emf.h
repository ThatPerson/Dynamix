#ifndef EMF_H
#define EMF_H

#undef MOD_EMF
#undef MOD_EMFT
#undef MOD_DEMF
#undef MOD_DEMFT
#define MOD_EMF 		2
#define MOD_EMFT 		3
#define MOD_DEMF		5 						///< Extended Model Free 
#define MOD_DEMFT		6	

inline Decimal J0_EMF(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f);
inline Decimal J0_EMF_CC(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f);
Decimal EMF_Dipolar_R1(Decimal omega_X, Decimal omega_Y, Decimal d, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f);
Decimal EMF_Dipolar_R2(Decimal omega_X, Decimal omega_Y, Decimal d, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal J0sum);
Decimal EMF_R1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, int mode);
Decimal EMF_R2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, int mode);

#include "emf.c"
#endif
