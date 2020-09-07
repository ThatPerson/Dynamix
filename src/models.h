#ifndef MODELS_H
#define MODELS_H

inline long double J0_SMF(long double omega, long double tau, long double S2);
inline long double J0_EMF(long double omega, long double taus, long double S2s, long double tauf, long double S2f);
inline long double J0_EMF_CC(long double omega, long double taus, long double S2s, long double tauf, long double S2f);
long double sq(long double x);
int sq_i(int x);
inline double GAF_Dipolar_R1(long double omega_obs, long double omega_neigh, long double taus, long double S2s, long double tauf, long double S2f, long double D);

inline double GAF_CSA_R2(long double omega, \
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
inline double GAF_Dipolar_R2(long double omega_obs, long double omega_neigh, double w1, double wr, long double taus, long double S2s, long double tauf, long double S2f, double D);
double SMF_R1(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, unsigned int mode);
double SMF_R2(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, unsigned int mode);
double EMF_R1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode);
double EMF_R2(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode);
int GAF_S2(long double sig[3], struct Orient ** A, struct Orient ** B, double * S2[], int length,  int mode);
double GAF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double GAF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double GAF_13CR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
double GAF_13CR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf);
#include "models.c"
#endif
