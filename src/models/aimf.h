#ifndef AIMF_H
#define AIMF_H

//#ifndef EMF_H
//#include "models/emf.h"
//#endif

#undef MOD_AIMF
#undef MOD_AIMFT
#define MOD_AIMF 			23
#define MOD_AIMFT 			24

double AIMF_Dipolar_R1(long double omega_obs, long double omega_neigh, long double taus, long double S2s, long double tauf, long double S2f, long double D);
double AIMF_CSA_R2(long double omega, \
						 double w1, \
						 double wr, \
						 long double taus, \
						 long double S2s, \
						 long double tauf, \
						 long double S2f, \
						 double D2, \
						 long double  (*J_SD)(\
						 	long double,\
						 	long double, \
						 	long double, \
						 	long double, \
						 	long double)\
						);
double AIMF_Dipolar_R2(long double omega_obs, long double omega_neigh, double w1, double wr, long double taus, long double S2s, long double tauf, long double S2f, double D);
int AIMF_S2(long double order_params[3], struct Orient ** A, double * S2[], int length);
double AIMF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double AIMF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double AIMF_13CR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double AIMF_13CR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);

#include "aimf.c"
#endif
