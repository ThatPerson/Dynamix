#ifndef SMF_H
#define SMF_H

#undef MOD_SMF
#undef MOD_SMFT
#define MOD_SMF 1
#define MOD_SMFT 4

inline double J0_SMF(double omega, double tau, double S2);
double SMF_Dipolar_R1(double omega_X, double omega_Y, double d, double tau, double S2);
double SMF_Dipolar_R2(double omega_X, double omega_Y, double d, double tau, double S2, double J0sum);

double SMF_R1(struct Residue *res, struct Relaxation* relax, double tau, double S2, unsigned int mode);
double SMF_R2(struct Residue *res, struct Relaxation* relax, double tau, double S2, unsigned int mode);

#include "smf.c"
#endif
