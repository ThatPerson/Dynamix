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
 * Calculates dipolar contribution to R1 relaxation rate.
 * @param omega_X
 *  Nuclei being observed
 * @param omega_Y
 *  Other nuclei (eg 1H)
 * @param d
 *  Dipolar coupling between nuclei
 * @param tau
 *  Correlation time of motion
 * @param S2
 *  Order parameter of motion.
 */
double SMF_Dipolar_R1(long double omega_X, long double omega_Y, long double d, long double tau, long double S2) {
	// omega_X is the measured nuclei (eg 15N, 13C), Y is the other (eg 1H)
	long double Jcomp = 0;
	Jcomp += J0_SMF(omega_Y - omega_X, tau, S2);
	Jcomp += 3. * J0_SMF(omega_X, tau, S2);
	Jcomp += 6. * J0_SMF(omega_Y + omega_X, tau, S2);
	return (double) (0.1 * d * d * Jcomp);
}


/**
 * Calculate R1p/R2 contribution from dipolar interaction
 * @param omega_X
 *  Observed nuclei
 * @param omega_Y
 *  Coupled nuclei
 * @param d
 *  Dipolar coupling
 * @param tau
 *  Characteristic time scale
 * @param S2
 *  Order parameter
 * @param J0sum
 *  J0sum
 */
double SMF_Dipolar_R2(long double omega_X, long double omega_Y, long double d, long double tau, long double S2, long double J0sum) {
	// all args should've been T_DOWN'ed
	long double JD = J0sum + 3 * J0_SMF(omega_X, tau, S2);
	JD += J0_SMF(omega_Y - omega_X, tau, S2);
	JD += 6 * J0_SMF(omega_Y, tau, S2);
	JD += 6 * J0_SMF(omega_Y + omega_X, tau, S2);
	return (double) ((1/20.) * d * d * JD);
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

	long double omega_1H = T_DOWN * 2 * M_PI * field;
	long double omega_15N = T_DOWN * 2 * M_PI * field / 9.869683408806043;
	long double omega_13C = T_DOWN * 2 * M_PI * field / 3.976489314034722;
	long double omega_L;
	long double d2tot;
	long double wCOCa = 120 * omega_13C * 0.000001;
	double *csa;
	long double R1D = 0, R1CSA = 0;
	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		//d = -D_NH;
		// double SMF_Dipolar_R1(long double omega_X, long double omega_Y, long double d, long double tau, long double S2) {
		R1D += SMF_Dipolar_R1(omega_15N, omega_1H, -D_NH, tau, S2); // N-H
		R1D += SMF_Dipolar_R1(omega_15N, omega_1H, -D_HNr, tau, S2); // N-Hr
		R1D += SMF_Dipolar_R1(omega_15N, omega_13C, -D_CN, tau, S2); // N-C
		R1D += SMF_Dipolar_R1(omega_15N, omega_13C, -D_CaN, tau, S2); // N-Ca		
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		R1D += SMF_Dipolar_R1(omega_13C, omega_1H, -D_CH, tau, S2); // C-H
		R1D += SMF_Dipolar_R1(omega_13C, omega_1H, -D_CHr, tau, S2); // C-Hr
		R1D += SMF_Dipolar_R1(omega_13C, omega_15N, -D_CN, tau, S2); // C-N
		R1D += SMF_Dipolar_R1(omega_13C, omega_13C + wCOCa, -D_CC, tau, S2); // C-C
		//d = -D_CH;
	} else {
		printf("Incorrect mode.\n");
		return -1;
	}

	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	omega_L *= T_DOWN;
	
	R1CSA = (2/15.) * d2tot * J0_SMF(omega_L, tau, S2);
	
	R1D = R1D * T_DOWN;
	R1CSA = R1CSA * T_DOWN;
	return (double) R1D + (double) R1CSA;
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

	long double omega_1H = T_DOWN * 2 * M_PI * field;
	long double omega_15N = T_DOWN * 2 * M_PI * field / 9.869683408806043;
	long double omega_13C = T_DOWN * 2 * M_PI * field / 3.976489314034722;
	long double omega_L;
	long double d2tot;
	double *csa;
	long double R2D = 0, R2CSA = 0;
	long double wCOCa = 120 * omega_13C * 0.000001;
	long double w1 = relax->w1;
	long double wr = relax->wr;
	w1 *= T_DOWN;
	wr *= T_DOWN;
	
	long double J0sum = 0;
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 - 2 * wr), tau, S2);
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 + 2 * wr), tau, S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 - wr), tau, S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 + wr), tau, S2);
	
	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		//d = -D_NH;
		R2D += SMF_Dipolar_R2(omega_15N, omega_1H, -D_NH, tau, S2, J0sum); // N-H
		R2D += SMF_Dipolar_R2(omega_15N, omega_1H, -D_HNr, tau, S2, J0sum); // N-Hr
		R2D += SMF_Dipolar_R2(omega_15N, omega_13C, -D_CN, tau, S2, J0sum); // N-C
		R2D += SMF_Dipolar_R2(omega_15N, omega_13C, -D_CaN, tau, S2, J0sum); // N-Ca	
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		R2D += SMF_Dipolar_R2(omega_13C, omega_1H, -D_CH, tau, S2, J0sum); // C-H
		R2D += SMF_Dipolar_R2(omega_13C, omega_1H, -D_CHr, tau, S2, J0sum); // C-Hr
		R2D += SMF_Dipolar_R2(omega_13C, omega_15N, -D_CN, tau, S2, J0sum); // C-N
		R2D += SMF_Dipolar_R2(omega_13C, omega_13C + wCOCa, -D_CC, tau, S2, J0sum); // C-C
		//d = -D_CH;
	} else {
		printf("Incorrect mode.\n");
		return -1;
	}
	

	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);
	
	omega_L *= T_DOWN;

	R2CSA = (1/45.) * d2tot * (J0sum + 3 * J0_SMF(omega_L, tau, S2));
	R2D *= T_DOWN;
	R2CSA *= T_DOWN;
	return (double) R2D + (double) R2CSA;
}
