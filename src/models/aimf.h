#ifndef AIMF_H
#define AIMF_H

//#ifndef EMF_H
//#include "models/emf.h"
//#endif

#undef MOD_AIMF
#undef MOD_AIMFT
#define MOD_AIMF 			23
#define MOD_AIMFT 			24

double AIMF_Dipolar_R1(double omega_obs, double omega_neigh, double taus, double S2s, double tauf, double S2f, double D);
double AIMF_CSA_R2(double omega, \
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
double AIMF_Dipolar_R2(double omega_obs, double omega_neigh, double w1, double wr, double taus, double S2s, double tauf, double S2f, double D);
int AIMF_S2(double order_params[3], struct Orient ** A, double * S2[], int length);
double AIMF_15NR1(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * sigs, double * sigf);
double AIMF_15NR2(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * sigs, double * sigf);
double AIMF_13CR1(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * sigs, double * sigf);
double AIMF_13CR2(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * sigs, double * sigf);

#include "aimf.c"
#endif
