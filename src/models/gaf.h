#ifndef GAF_H
#define GAF_H

#ifndef EMF_H
#include "models/emf.h"
#endif

#undef MOD_GAF
#undef MOD_GAFT
#define MOD_GAF 			7
#define MOD_GAFT 			8

double GAF_Dipolar_R1(long double omega_obs, long double omega_neigh, long double taus, long double S2s, long double tauf, long double S2f, long double D);
double GAF_CSA_R2(long double omega, \
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
double GAF_Dipolar_R2(long double omega_obs, long double omega_neigh, double w1, double wr, long double taus, long double S2s, long double tauf, long double S2f, double D);
int GAF_S2(long double sig[3], struct Orient ** A, struct Orient ** B, double * S2[], int length,  int mode);
double GAF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double GAF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double GAF_13CR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double GAF_13CR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);

#include "gaf.c"
#endif
