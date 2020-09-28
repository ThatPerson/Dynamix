#ifndef SMF_H
#define SMF_H

#undef MOD_SMF
#undef MOD_SMFT
#define MOD_SMF 1
#define MOD_SMFT 4

inline long double J0_SMF(long double omega, long double tau, long double S2);
double SMF_Dipolar_R1(long double omega_X, long double omega_Y, long double d, long double tau, long double S2);
double SMF_Dipolar_R2(long double omega_X, long double omega_Y, long double d, long double tau, long double S2, long double J0sum);

double SMF_R1(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, unsigned int mode);
double SMF_R2(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, unsigned int mode);

#include "smf.c"
#endif
