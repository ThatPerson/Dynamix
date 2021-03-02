/**
 * @file egaf.c
 */

/* There are undoubtedly much better ways of doing this that don't just involve copying the GAF code.
 * I may come back to it if I have time.
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>


/**
 * Calculates R1 relaxation rate for 15N nucleus under gaussian axial fluctuations for slow motions, and order parameter for fast motions.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @warning I am unsure of the form of the spectral density for cross correlated (eg Jxy) motions.
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param sigs
 *  Deflection angles [alpha, beta, gamma] for slow motion
 * @param S2f
 *  Order parameter for fast motion
 * @return R1
 *  Returns R1 as double
 */
double EGAF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double S2f) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N;
	double d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;

	double *csa;
	csa = res->csaN;

	d2x = (double) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = (double) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= (double) sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* Calculate all order parameters */
	double S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	double *S2s[] = {&S2NHs, &S2NCSAxs, &S2NCSAys, &S2NCSAxys, &S2CNs, &S2CaNs};
	
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	struct Orient * Bs[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);

	/* N CSA relaxation contribution */
	long double R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1NH, R1NHr, R1CN, R1CaN;
	long double J1 = 0;
	

	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	

	J1 = J0_EMF(omega_15N, taus, S2NCSAxs, tauf, S2f);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
	J1 = J0_EMF(omega_15N, taus, S2NCSAys, tauf, S2f);
	R1CSAy = (1/15.) * d2y * J1;
	J1 = J0_EMF_CC(omega_15N, taus, S2NCSAxys, tauf, S2f);
	R1CSAxy = (1/15.) * d2xy * J1;
	
	/** Eq A30, Bremi1997 */
	R1CSA = R1CSAx + R1CSAy + 2.*R1CSAxy;
	



	/* N Dipolar Interactions Contributions */
	R1NH = GAF_Dipolar_R1(omega_15N, omega_1H, taus, S2NHs, tauf, S2f, D_NH);
	R1NHr = GAF_Dipolar_R1(omega_15N, omega_1H, taus, S2CNs, tauf, S2f, D_HNr);
	R1CN = GAF_Dipolar_R1(omega_15N, omega_13C, taus, S2CNs, tauf, S2f, D_CN);
	R1CaN = GAF_Dipolar_R1(omega_15N, omega_13C, taus, S2CaNs, tauf, S2f, D_CaN);

	return (double) (R1CSA + R1NH + R1NHr + R1CN + R1CaN) * T_DOWN;
}

/**
 * Calculates R2 relaxation rate for 15N nucleus under (slow) gaussian axial fluctuations and isotropic fast motions.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param sigs
 *  Deflection angles [alpha, beta, gamma] for slow motion
 * @param S2f
 *  order parameter for fast motion
 * @return R2
 *  Returns R2 as double
 */
double EGAF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double S2f) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz
	double w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N;
	double d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;

	double *csa;
	csa = res->csaN;

	d2x = (double) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = (double) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= (double) sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	double S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	
	double *S2s[] = {&S2NHs, &S2NCSAxs, &S2NCSAys, &S2NCSAxys, &S2CNs, &S2CaNs};
	
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	struct Orient * Bs[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);

	/* N CSA relaxation contribution */
	double R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2NH, R2NHr, R2CN, R2CaN;

	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	R2CSAx = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAxs, tauf, S2f, d2x, J0_EMF);
	R2CSAy = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAys, tauf, S2f, d2y, J0_EMF);
	R2CSAxy = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAxys, tauf, S2f, d2xy, J0_EMF_CC);
	/**
	 * Equation for spectral density for xy taken from Manuscript
	 */ 
	R2CSA = R2CSAx + R2CSAy + 2.* R2CSAxy;

	/* N Dipolar Interactions Contributions */
	R2NH = GAF_Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, S2NHs, tauf, S2f, D_NH);
	R2NHr = GAF_Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, S2CNs, tauf, S2f, D_HNr);
	R2CN = GAF_Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, S2CNs, tauf, S2f, D_CN);
	R2CaN = GAF_Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, S2CaNs, tauf, S2f, D_CaN);
	return (R2CSA + R2NH + R2NHr + R2CN + R2CaN)*T_DOWN;
}



 
 /**
 * Calculates R1 relaxation rate for 13C nucleus under gaussian axial fluctuations.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param sigs
 *  Deflection angles [alpha, beta, gamma] for slow motion
 * @param sigf
 *  Deflection angles [alpha, beta, gamma] for fast motion
 * @return R1
 *  Returns R1 as double
 */

double EGAF_13CR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double S2f) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N, wCOCa;
	double d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	double *csa;
	csa = res->csaC;

	d2x = (double) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
	d2y = (double) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
	d2xy= (double) sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* Calculate all order parameters */
	double S2CHs, S2CSAxs, S2CSAys, S2CSAxys, S2CNs, S2CCs;
	
	double *S2s[] = {&S2CHs, &S2CSAxs, &S2CSAys, &S2CSAxys, &S2CNs, &S2CCs};
	
	struct Orient * As[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	struct Orient * Bs[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	
	/* N CSA relaxation contribution */
	long double R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1CH, R1CHr, R1CN, R1CC;
	long double J1 = 0;


	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	J1 = J0_EMF(omega_13C, taus, S2CSAxs, tauf, S2f);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
	J1 = J0_EMF(omega_13C, taus, S2CSAys, tauf, S2f);
	R1CSAy = (1/15.) * d2y * J1;
	J1 = J0_EMF_CC(omega_13C, taus, S2CSAxys, tauf, S2f);
	R1CSAxy = (1/15.) * d2xy * J1;
	R1CSA = R1CSAx + R1CSAy + 2.*R1CSAxy;
	
	

	/* N Dipolar Interactions Contributions */
	R1CH = GAF_Dipolar_R1(omega_13C, omega_1H, taus, S2CHs, tauf, S2f, D_CH);
	R1CHr = GAF_Dipolar_R1(omega_13C, omega_1H, taus, S2CNs, tauf, S2f, D_CHr);
	R1CN = GAF_Dipolar_R1(omega_13C, omega_15N, taus, S2CNs, tauf, S2f, D_CN);
	R1CC = GAF_Dipolar_R1(omega_13C, omega_13C - wCOCa, taus, S2CCs, tauf, S2f, D_CC);

	return (double) (R1CSA + R1CH + R1CHr + R1CN + R1CC)*T_DOWN;
}

/**
 * Calculates R2 relaxation rate for 13C nucleus under (slow) gaussian axial fluctuations and (fast) isotropic motions.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param sigs
 *  Deflection angles [alpha, beta, gamma] for slow motion
 * @param S2f
 *  Order parameter for fast motion
 * @return R2
 *  Returns R2 as double
 */

double EGAF_13CR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double S2f) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz
	double w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N, wCOCa;
	double d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	double *csa;
	csa = res->csaC;

	d2x = (double) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
	d2y = (double) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
	d2xy= (double) sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	double S2CHs, S2CSAxs, S2CSAys, S2CSAxys, S2CNs, S2CCs;

	double *S2s[] = {&S2CHs, &S2CSAxs, &S2CSAys, &S2CSAxys, &S2CNs, &S2CCs};
	
	struct Orient * As[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	struct Orient * Bs[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	

	/* CSA relaxation contribution */
	double R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2CH, R2CHr, R2CN, R2CC;
	
	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	R2CSAx = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAxs, tauf, S2f, d2x, J0_EMF);
	R2CSAy = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAys, tauf, S2f, d2y, J0_EMF);
	R2CSAxy = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAxys, tauf, S2f, d2xy, J0_EMF_CC);
	R2CSA = R2CSAx + R2CSAy + 2.*R2CSAxy;

	/* N Dipolar Interactions Contributions */
	R2CH = GAF_Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, S2CHs, tauf, S2f, D_CH);
	R2CHr = GAF_Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, S2CNs, tauf, S2f, D_CHr);
	R2CN = GAF_Dipolar_R2(omega_13C, omega_15N, w1, wr, taus, S2CNs, tauf, S2f, D_CN);
	R2CC = GAF_Dipolar_R2(omega_13C, omega_13C - wCOCa, w1, wr, taus, S2CCs, tauf, S2f, D_CC);
	return (double) ((R2CSA + R2CH + R2CHr + R2CN + R2CC)*(long double)T_DOWN);
}
