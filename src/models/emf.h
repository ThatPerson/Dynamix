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

inline long double J0_EMF(long double omega, long double taus, long double S2s, long double tauf, long double S2f);
inline long double J0_EMF_CC(long double omega, long double taus, long double S2s, long double tauf, long double S2f);
double EMF_Dipolar_R1(long double omega_X, long double omega_Y, long double d, long double taus, long double S2s, long double tauf, long double S2f);
double EMF_Dipolar_R2(long double omega_X, long double omega_Y, long double d, long double taus, long double S2s, long double tauf, long double S2f, long double J0sum);
double EMF_R1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode);
double EMF_R2(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode);

#include "emf.c"
#endif
