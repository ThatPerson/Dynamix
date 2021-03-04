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
inline double J0_EMF(double omega, double taus, double S2s, double tauf, double S2f) {
	return ( \
			((((1 - (double) S2f)) * (double) tauf)) \
			/ (1 + ((double) omega * (double) tauf * (double) omega * (double) tauf)) \
		)\
		+\
		(\
			(((double) S2f) * (1 - (double) S2s) * (double) taus)\
			/ (1 + ((double) omega * (double) taus * (double) omega * (double) taus))\
		);
}

/** Cross correlated spectral density. Taken from Manuscript, 
 * Juv(w) = (1-S2f) tf / (1 + (wtf)^2) + (1/P2) S2f (P2 - S2s) (ts / 1+(wts)^2)
 * eq 2
 */
inline double J0_EMF_CC(double omega, double taus, double S2s, double tauf, double S2f) {
	return ( \
			((((1 - (double) S2f)) * (double) tauf)) \
			/ (1 + ((double) omega * (double) tauf * (double) omega * (double) tauf)) \
		)\
		-\
		2. * (\
			(((double) S2f) * (-0.5 - (double) S2s) * (double) taus)\
			/ (1 + ((double) omega * (double) taus * (double) omega * (double) taus))\
		);
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
double EMF_Dipolar_R1(double omega_X, double omega_Y, double d, double taus, double S2s, double tauf, double S2f) {
	// omega_X is the measured nuclei (eg 15N, 13C), Y is the other (eg 1H)
	double Jcomp = 0;
	Jcomp += J0_EMF(omega_Y - omega_X, taus, S2s, tauf, S2f);
	Jcomp += 3. * J0_EMF(omega_X, taus, S2s, tauf, S2f);
	Jcomp += 6. * J0_EMF(omega_Y + omega_X, taus, S2s, tauf, S2f);
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
double EMF_Dipolar_R2(double omega_X, double omega_Y, double d, double taus, double S2s, double tauf, double S2f, double J0sum) {
	// all args should've been T_DOWN'ed
	double JD = J0sum + 3 * J0_EMF(omega_X, taus, S2s, tauf, S2f);
	JD += J0_EMF(omega_Y - omega_X, taus, S2s, tauf, S2f);
	JD += 6 * J0_EMF(omega_Y, taus, S2s, tauf, S2f);
	JD += 6 * J0_EMF(omega_Y + omega_X, taus, S2s, tauf, S2f);
	return (double) ((1/20.) * d * d * JD);
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
double EMF_R1(struct Residue *res, struct Relaxation* relax, double taus, double S2s, double tauf, double S2f, int mode) {
	double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a double128, and
	 * software division) I've predefined it. */

	double omega_1H = T_DOWN * 2 * M_PI * field;
	double omega_15N = T_DOWN * 2 * M_PI * field / 9.869683408806043;
	double omega_13C = T_DOWN * 2 * M_PI * field / 3.976489314034722;
	double wCOCa = 120 * omega_13C * 0.000001;
	double omega_L;
	double d2tot;
	double *csa;
	
	double R1D = 0, R1CSA = 0;

	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		//d = -D_NH;
		R1D += EMF_Dipolar_R1(omega_15N, omega_1H, -D_NH, taus, S2s, tauf, S2f); // N-H
		R1D += EMF_Dipolar_R1(omega_15N, omega_1H, -D_HNr, taus, S2s, tauf, S2f); // N-Hr
		R1D += EMF_Dipolar_R1(omega_15N, omega_13C, -D_CN, taus, S2s, tauf, S2f); // N-C
		R1D += EMF_Dipolar_R1(omega_15N, omega_13C, -D_CaN, taus, S2s, tauf, S2f); // N-Ca	
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		//d = -D_CH;
		R1D += EMF_Dipolar_R1(omega_13C, omega_1H, -D_CH, taus, S2s, tauf, S2f); // C-H
		R1D += EMF_Dipolar_R1(omega_13C, omega_1H, -D_CHr, taus, S2s, tauf, S2f); // C-Hr
		R1D += EMF_Dipolar_R1(omega_13C, omega_15N, -D_CN, taus, S2s, tauf, S2f); // C-N
		R1D += EMF_Dipolar_R1(omega_13C, omega_13C - wCOCa, -D_CC, taus, S2s, tauf, S2f); // C-C
	} else {
		ERROR("Incorrect mode.");
		return -1;
	}

	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	omega_L *= T_DOWN;

	R1CSA = (2/15.) * d2tot * J0_EMF(omega_L, taus, S2s, tauf, S2f);

	R1D *= T_DOWN;
	R1CSA *= T_DOWN;
	return (double) R1D + (double) R1CSA;
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
double EMF_R2(struct Residue *res, struct Relaxation* relax, double taus, double S2s, double tauf, double S2f, int mode) {
	//if (d_mode == DIPOLAR) 
	//	S2f = res->S2_dipolar / S2s;
	double field = relax->field * 1000000; // conversion to Hz
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a double128, and
	 * software division) I've predefined it.  */

	double omega_1H = T_DOWN * 2 * M_PI * field;
	double omega_15N = T_DOWN * 2 * M_PI * field / 9.869683408806043;
	double omega_13C = T_DOWN * 2 * M_PI * field / 3.976489314034722;
	double omega_L;
	double d2tot;
	double *csa;
	double wCOCa = 120 * omega_13C * 0.000001;
	
	double w1 = relax->w1;
	double wr = relax->wr;
	w1 *= T_DOWN;
	wr *= T_DOWN;
	
	double J0sum = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f);

	double R2D = 0, R2CSA = 0;
	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		//d = -D_NH;
		R2D += EMF_Dipolar_R2(omega_15N, omega_1H, -D_NH, taus, S2s, tauf, S2f, J0sum); // N-H
		R2D += EMF_Dipolar_R2(omega_15N, omega_1H, -D_HNr, taus, S2s, tauf, S2f, J0sum); // N-Hr
		R2D += EMF_Dipolar_R2(omega_15N, omega_13C, -D_CN, taus, S2s, tauf, S2f, J0sum); // N-C
		R2D += EMF_Dipolar_R2(omega_15N, omega_13C, -D_CaN, taus, S2s, tauf, S2f, J0sum); // N-Ca	
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		//d = -D_CH;
		R2D += EMF_Dipolar_R2(omega_13C, omega_1H, -D_CH, taus, S2s, tauf, S2f, J0sum); // C-H
		R2D += EMF_Dipolar_R2(omega_13C, omega_1H, -D_CHr, taus, S2s, tauf, S2f, J0sum); // C-Hr
		R2D += EMF_Dipolar_R2(omega_13C, omega_15N, -D_CN, taus, S2s, tauf, S2f, J0sum); // C-N
		R2D += EMF_Dipolar_R2(omega_13C, omega_13C - wCOCa, -D_CC, taus, S2s, tauf, S2f, J0sum); // C-C
	} else {
		ERROR("Incorrect mode.");
		return -1;
	}

	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	omega_L *= T_DOWN;

	R2D *= T_DOWN;
	R2CSA = (1/45.) * d2tot * (J0sum + 3 * J0_EMF(omega_L, taus, S2s, tauf, S2f)) * T_DOWN;

	return (double) R2D + (double) R2CSA;
}

