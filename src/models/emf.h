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

inline double J0_EMF(double omega, double taus, double S2s, double tauf, double S2f);
inline double J0_EMF_CC(double omega, double taus, double S2s, double tauf, double S2f);
double EMF_Dipolar_R1(double omega_X, double omega_Y, double d, double taus, double S2s, double tauf, double S2f);
double EMF_Dipolar_R2(double omega_X, double omega_Y, double d, double taus, double S2s, double tauf, double S2f, double J0sum);
double EMF_R1(struct Residue *res, struct Relaxation* relax, double taus, double S2s, double tauf, double S2f, int mode);
double EMF_R2(struct Residue *res, struct Relaxation* relax, double taus, double S2s, double tauf, double S2f, int mode);

#include "emf.c"
#endif
