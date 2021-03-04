#ifndef GAF_H
#define GAF_H

#ifndef EMF_H
#include "models/emf.h"
#endif

#undef MOD_GAF
#undef MOD_GAFT
#define MOD_GAF 			7
#define MOD_GAFT 			8

double GAF_Dipolar_R1(double omega_obs, double omega_neigh, double taus, double S2s, double tauf, double S2f, double D);
double GAF_CSA_R2(double omega, \
						 double w1, \
						 double wr, \
						 double taus, \
						 double S2s, \
						 double tauf, \
						 double S2f, \
						 double D2, \
						 double  (*J_SD)(\
						 	double,\
						 	double, \
						 	double, \
						 	double, \
						 	double)\
						);
double GAF_Dipolar_R2(double omega_obs, double omega_neigh, double w1, double wr, double taus, double S2s, double tauf, double S2f, double D);
int GAF_S2(double sig[3], struct Orient ** A, struct Orient ** B, double * S2[], int length,  int mode);
double GAF_15NR1(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * sigs, double * sigf);
double GAF_15NR2(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * sigs, double * sigf);
double GAF_13CR1(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * sigs, double * sigf);
double GAF_13CR2(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * sigs, double * sigf);

#include "gaf.c"
#endif
