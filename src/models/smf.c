/**
 * @file smf.c
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>


/**
 * Calculates spectral density function for given frequency (omega) according to simple model free analysis.
 * Implements J(w) = (1 - S2) tau / (1 + (w tau)^2) as from Lipari1982 eq 35 (assuming no overall tumbling)
 * and the factor of 2/5 has been taken into the summation (see below).
 */
inline long double J0_SMF(long double omega, long double tau, long double S2) {
	return (((1 - S2) * tau)) \
		/ (1 + (omega * omega * tau * tau));
}


/**
 * Calculates R1 relaxation rate using simple model free analysis.
 * Kurbanov2011, eqs 8-10
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param tau
 *  Correlation time of motion being optimized
 * @param S2
 *  Order parameter of motion
 * @param mode
 *  One of MODE_15N or MODE_13C depending on relaxation data type considered.
 * @return R1
 *  Returns R1 as double
 */
double SMF_R1(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, unsigned int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */

	long double field = relax->field * 1000000; // conversion to Hz
	//tau = tau * T_UP;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Constants are in datatypes.c */

	long double omega_1H = 2 * M_PI * field;
	long double omega_L;
	long double d = 0, d2tot;
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
	
	omega_L *= T_DOWN;
	omega_1H *= T_DOWN;
	double R1D = 0, R1CSA = 0;
	long double Jcomp = 0;
	Jcomp += J0_SMF(omega_1H - omega_L, tau, S2);
	Jcomp += 3 * J0_SMF(omega_L, tau, S2);
	Jcomp += 6 * J0_SMF(omega_1H + omega_L, tau, S2);
	R1D = 0.1 * d * d * Jcomp;
	R1CSA = (2/15.) * d2tot * J0_SMF(omega_L, tau, S2);
	R1D = R1D * T_DOWN;
	R1CSA = R1CSA * T_DOWN;
	return R1D + R1CSA;
}


/**
 * Calculates R2 relaxation rate using simple model free analysis
 * Kurbanov2011, eqs 8-10
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param tau
 *  Correlation time of motion being optimized
 * @param S2
 *  Order parameter of motion
 * @param mode
 *  One of MODE_15N or MODE_13C depending on relaxation data type considered.
 * @return R2
 *  Returns R2 as double
 */
double SMF_R2(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, unsigned int mode) {

	long double field = relax->field * 1000000; // conversion to Hz

	//tau *= T_UP;
	
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
	
	long double w1 = relax->w1;
	long double wr = relax->wr;

	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_L *= T_DOWN;
	
	long double R2NH = 0, R2NCSA = 0;
	long double J0sum = 0;
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 - 2 * wr), tau, S2);
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 + 2 * wr), tau, S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 - wr), tau, S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 + wr), tau, S2);

	long double JNH = J0sum + 3 * J0_SMF(omega_L, tau, S2);
	JNH += J0_SMF(omega_1H - omega_L, tau, S2);
	JNH += 6 * J0_SMF(omega_1H, tau, S2);
	JNH += 6 * J0_SMF(omega_1H + omega_L, tau, S2);


	R2NH = (1/20.) * d * d * JNH;
	R2NCSA = (1/45.) * d2tot * (J0sum + 3 * J0_SMF(omega_L, tau, S2));
	R2NH *= T_DOWN;
	R2NCSA *= T_DOWN;
	return R2NH + R2NCSA;
}

