/**
 * @file emf.c
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
  

/** Calculates spectral density function for given frequency according to extended model free analysis.
 *  Implements J(w) = (1 - S2f) tauf / (1 + (w tauf)^2) + S2f (1 - S2s) taus / (1 + (w taus)^2), 
 *  as in eq 2, Clore 1990 (assuming S2 = S2s*S2f and that there is no overall tumbling).
 */
inline long double J0_EMF(long double omega, long double taus, long double S2s, long double tauf, long double S2f) {
	return ( \
			((((1 - (double) S2f)) * (long double) tauf)) \
			/ (1 + ((double) omega * (long double) tauf * (double) omega * (long double) tauf)) \
		)\
		+\
		(\
			(((double) S2f) * (1 - (double) S2s) * (long double) taus)\
			/ (1 + ((double) omega * (long double) taus * (double) omega * (long double) taus))\
		);
}

/** Cross correlated spectral density. Taken from Manuscript, 
 * Juv(w) = (1-S2f) tf / (1 + (wtf)^2) + (1/P2) S2f (P2 - S2s) (ts / 1+(wts)^2)
 * eq 2
 */
inline long double J0_EMF_CC(long double omega, long double taus, long double S2s, long double tauf, long double S2f) {
	return ( \
			((((1 - (double) S2f)) * (long double) tauf)) \
			/ (1 + ((double) omega * (long double) tauf * (double) omega * (long double) tauf)) \
		)\
		-\
		2. * (\
			(((double) S2f) * (-0.5 - (double) S2s) * (long double) taus)\
			/ (1 + ((double) omega * (long double) taus * (double) omega * (long double) taus))\
		);
}
  

/**
 * Calculates R1 relaxation rate using extended model free analysis
 * Kurbanov2011, 8-10
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param S2s
 *  Order parameter of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param S2f
 *  Order parameter of fast motion
 * @param mode
 *  One of MODE_15N or MODE_13C depending on relaxation data type considered.
 * @return R1
 *  Returns R1 as double
 */
double EMF_R1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode) {
	long double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. */

	long double omega_1H = 2 * M_PI * field;
	long double omega_L;
	long double d, d2tot;
	double *csa;

	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		d = -D_NH;
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		d = -D_CH;
	} else {
		printf("Incorrect mode.\n");
		return -1;
	}

	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	double R1D = 0, R1CSA = 0;
	long double Jcomp = 0;
	
	omega_1H *= T_DOWN;
	omega_L *= T_DOWN;

	Jcomp += J0_EMF(omega_1H - omega_L, taus, S2s, tauf, S2f);
	Jcomp += 3 * J0_EMF(omega_L, taus, S2s, tauf, S2f);
	Jcomp += 6 * J0_EMF(omega_1H + omega_L, taus, S2s, tauf, S2f);
	R1D = 0.1 * d * d * Jcomp;
	R1CSA = (2/15.) * d2tot * J0_EMF(omega_L, taus, S2s, tauf, S2f);

	R1D *= T_DOWN;
	R1CSA *= T_DOWN;
	return R1D + R1CSA;
}

/**
 * Calculates R2 relaxation rate using extended model free analysis
 * Kurbanov2011, 8-10
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param S2s
 *  Order parameter of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param S2f
 *  Order parameter of fast motion.
 * @param mode
 *  One of MODE_15N or MODE_13C depending on relaxation data type considered.
 * @return R2
 *  Returns R2 as double
 */
double EMF_R2(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode) {
	//if (d_mode == DIPOLAR) 
	//	S2f = res->S2_dipolar / S2s;
	long double field = relax->field * 1000000; // conversion to Hz
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it.  */

	long double omega_1H = 2 * M_PI * field;
	long double omega_L;
	long double d, d2tot;
	double *csa;

	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		d = -D_NH;
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		d = -D_CH;
	} else {
		printf("Incorrect mode.\n");
		return -1;
	}

	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	long double w1 = relax->w1;
	long double wr = relax->wr;
	
	omega_1H *= T_DOWN;
	omega_L *= T_DOWN;
	w1 *= T_DOWN;
	wr *= T_DOWN;

	long double R2D = 0, R2CSA = 0;
	long double J0sum = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f);

	long double JNH = J0sum + 3 * J0_EMF(omega_L, taus, S2s, tauf, S2f);
	JNH += J0_EMF(omega_1H - omega_L, taus, S2s, tauf, S2f);
	JNH += 6 * J0_EMF(omega_1H, taus, S2s, tauf, S2f);
	JNH += 6 * J0_EMF(omega_1H + omega_L, taus, S2s, tauf, S2f);


	R2D = (1/20.) * d * d * JNH * T_DOWN;
	R2CSA = (1/45.) * d2tot * (J0sum + 3 * J0_EMF(omega_L, taus, S2s, tauf, S2f)) * T_DOWN;

	return R2D + R2CSA;
}

