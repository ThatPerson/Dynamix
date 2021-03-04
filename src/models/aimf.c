/**
 * @file aimf.c
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>



  
/** 
 * Calculates R1 contribution under AIMF model for Dipolar interaction between two atoms.
 * omega_obs is the frequency of the observed nucleus, omega_neigh is the frequency of the other atoms.\n
 * D is the dipolar coupling in rad/s.
 * Equation as 3.28, McConnell1986, only for solid state spins.
 * Kurbanov2011, eq 9 (factor of 1/10 arising from (2/5.) in spectral density and 1/4. in dipolar.
 */
double AIMF_Dipolar_R1(double omega_obs, double omega_neigh, double taus, double S2s, double tauf, double S2f, double D) {
	double q = (0.1) * sq(D) * ( \
		(J0_EMF(omega_neigh - omega_obs, taus, S2s, tauf, S2f)) \
		 + 3 * (J0_EMF(omega_obs, taus, S2s, tauf, S2f)) \
		  + 6 * (J0_EMF(omega_neigh + omega_obs, taus, S2s, tauf, S2f)) \
	  );
	return (double) q;
}
    

/**
 * Calculates the R2 contribution of CSA of nuclei with frequency omega.
 * D2 is the squared D22x/D22y/D22xy variable (dependent on CSA)
 * Kurbanov2011, eqs 18-20.
 */
 //double J0_EMF_CC(double omega, double taus, double S2s, double tauf, double S2f) {
 
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
						) {
	return (double) ( \
  		(1/45.) * (double) D2 * ( \
  			(2/3.) * J_SD(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f) +\
  			(2/3.) * J_SD(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J_SD(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J_SD(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f) + \
  			(3.)   * J_SD(omega, taus, S2s, tauf, S2f)\
  		)\
  	);
}

/** 
 * Calculates R2 contribution of dipolar interaction between two atoms as in AIMF_Dipolar_R1
 * Kurbanov2011, for the case where theta_p = 90 such that sin(theta_p) = 1 and so
 * R1p = 0.5*R1 + R1delta
 */
double AIMF_Dipolar_R2(double omega_obs, double omega_neigh, double w1, double wr, double taus, double S2s, double tauf, double S2f, double D) {
	return (double) (\
  		(1/20.) * sq(D) * (\
  			(2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f) +\
  			(2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f) +\
  			(3.)   * J0_EMF(omega_obs, taus, S2s, tauf, S2f) +\
  			(1.)   * J0_EMF(omega_neigh - omega_obs, taus, S2s, tauf, S2f) +\
  			(6.)   * J0_EMF(omega_neigh, taus, S2s, tauf, S2f) +\
  			(6.)   * J0_EMF(omega_neigh + omega_obs, taus, S2s, tauf, S2f)\
  		)\
  	);
}


/**
 * Calculates the order parameter for anisotropic model free motion. Essentially, the order_params define an ellipsoid.
 * The order parameter for a given theta, phi is the radius of the ellipsoid in that direction. eg;
 *   Let Sa, Sb, Sc be the squared order parameters along the following directions;
 *     Sa - CaCa axis (eg the axis with gamma deflections about it)
 *     Sb - ~N-H axis (eg the axis with alpha deflections about it)
 *     Sc - perpendicular to plane axis (eg beta).
 *   Then the equation for an ellipsoid is;
 *    (sin(theta)^2 sin(phi)^2 / Sa^2 + sin(theta)^2 cos(phi)^2 / Sb^2 + cos(theta)^2 / Sc^2 = 1 / r^2
 *   Evaluating this for r for a given theta, phi gives us the amplitude of order parameter in that orientation.
 * @param order_params
 *  Array containing [Sa, Sb, Sc] squared order parameters
 * @param A
 *  Array containing pointers to orientation vectors for 'A' for each S2 being calculated.
 * @param S2
 *  Array containing pointers to S2 values - this is where output will be put.
 * @param length
 *  Number of S2 values being calculated (eg length of S2)
 * @return Returns 1 in all cases.
 */
int AIMF_S2(double order_params[3], struct Orient ** A, double * S2[], int length) {
	int i;
	double ct, st, cp, sp;
	double S2a = powl(order_params[0], 2.), S2b = powl(order_params[1], 2.), S2c = powl(order_params[2], 2.);
	double SaSbSc = S2a * S2b * S2c;
	double denom;
	for (i = 0; i < length; i++) {
		// A[i].theta, A[i].phi
		// order_params[0], [1], [2] (x, y, z)
		ct = cos(A[i]->theta);
		st = sin(A[i]->theta);
		cp = cos(A[i]->phi);
		sp = sin(A[i]->phi);
		
		denom = 0;
		denom += S2b * S2c * powl(st * sp, 2.);
		denom += S2a * S2c * powl(st * cp, 2.);
		denom += S2a * S2b * powl(ct, 2.);
		
		*S2[i] = sqrtl(SaSbSc / denom);
		
	}
	return 1;
}





/**
 * Calculates R1 relaxation rate for 15N nucleus under anisotropic model free.
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param Ss
 *  Order parameters [alpha, beta, gamma] for slow motion
 * @param Sf
 *  Order parameters [alpha, beta, gamma] for fast motion
 * @return R1
 *  Returns R1 as double
 */
double AIMF_15NR1(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * Ss, double * Sf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a double128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N;
	double d2x, d2y;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;

	double *csa;
	csa = res->csaN;

	d2x = (double) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = (double) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	//d2xy= (double) sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* Calculate all order parameters */
	double S2NHs, S2NCSAxs, S2NCSAys, S2CNs, S2CaNs;
	double S2NHf, S2NCSAxf, S2NCSAyf, S2CNf, S2CaNf;
	
	double *S2s[] = {&S2NHs, &S2NCSAxs, &S2NCSAys, &S2CNs, &S2CaNs};
	double *S2f[] = {&S2NHf, &S2NCSAxf, &S2NCSAyf, &S2CNf, &S2CaNf};
	
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	
	AIMF_S2(Ss, As, S2s, 5);
	AIMF_S2(Sf, As, S2f, 5);
	

	/* N CSA relaxation contribution */
	double R1CSAx, R1CSAy, R1CSA, R1NH, R1NHr, R1CN, R1CaN;
	double J1 = 0;
	

	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	

	J1 = J0_EMF(omega_15N, taus, S2NCSAxs, tauf, S2NCSAxf);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
	J1 = J0_EMF(omega_15N, taus, S2NCSAys, tauf, S2NCSAyf);
	R1CSAy = (1/15.) * d2y * J1;
	//J1 = J0_EMF_CC(omega_15N, taus, S2NCSAxys, tauf, S2NCSAxyf);
	//R1CSAxy = (1/15.) * d2xy * J1;
	
	/** Eq A30, Bremi1997 */
	R1CSA = R1CSAx + R1CSAy;// + 2.*R1CSAxy;
	



	/* N Dipolar Interactions Contributions */
	R1NH = AIMF_Dipolar_R1(omega_15N, omega_1H, taus, S2NHs, tauf, S2NHf, D_NH);
	R1NHr = AIMF_Dipolar_R1(omega_15N, omega_1H, taus, S2CNs, tauf, S2CNf, D_HNr);
	R1CN = AIMF_Dipolar_R1(omega_15N, omega_13C, taus, S2CNs, tauf, S2CNf, D_CN);
	R1CaN = AIMF_Dipolar_R1(omega_15N, omega_13C, taus, S2CaNs, tauf, S2CaNf, D_CaN);

	return (double) (R1CSA + R1NH + R1NHr + R1CN + R1CaN) * T_DOWN;
}

/**
 * Calculates R2 relaxation rate for 15N nucleus under anisotropic model free.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param Ss
 *  Order parameters [a, b, c] for slow motion
 * @param Sf
 *  Order parameters [a, b, c] for fast motion
 * @return R2
 *  Returns R2 as double
 */
double AIMF_15NR2(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * Ss, double * Sf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz
	double w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a double128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N;
	double d2x, d2y;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;

	double *csa;
	csa = res->csaN;

	d2x = (double) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = (double) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);


	/* Calculate all order parameters */
	double S2NHs, S2NCSAxs, S2NCSAys, S2CNs, S2CaNs;
	double S2NHf, S2NCSAxf, S2NCSAyf, S2CNf, S2CaNf;
	
	double *S2s[] = {&S2NHs, &S2NCSAxs, &S2NCSAys, &S2CNs, &S2CaNs};
	double *S2f[] = {&S2NHf, &S2NCSAxf, &S2NCSAyf, &S2CNf, &S2CaNf};
	
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	
	AIMF_S2(Ss, As, S2s, 5);
	AIMF_S2(Sf, As, S2f, 5);
	
	

	/* N CSA relaxation contribution */
	double R2CSAx, R2CSAy, R2CSA, R2NH, R2NHr, R2CN, R2CaN;

	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	R2CSAx = AIMF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAxs, tauf, S2NCSAxf, d2x, J0_EMF);
	R2CSAy = AIMF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAys, tauf, S2NCSAyf, d2y, J0_EMF);
	//R2CSAxy = AIMF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAxys, tauf, S2NCSAxyf, d2xy, J0_EMF_CC);
	/**
	 * Equation for spectral density for xy taken from Manuscript
	 */ 
	R2CSA = R2CSAx + R2CSAy;// + 2.* R2CSAxy;

	/* N Dipolar Interactions Contributions */
	R2NH = AIMF_Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, S2NHs, tauf, S2NHf, D_NH);
	R2NHr = AIMF_Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, S2CNs, tauf, S2CNf, D_HNr);
	R2CN = AIMF_Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, S2CNs, tauf, S2CNf, D_CN);
	R2CaN = AIMF_Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, S2CaNs, tauf, S2CaNf, D_CaN);
	return (R2CSA + R2NH + R2NHr + R2CN + R2CaN)*T_DOWN;
}



 
 /**
 * Calculates R1 relaxation rate for 13C nucleus under anisotropic model free.

 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param Ss
 *  Order parameters [a, b, c] for slow motion
 * @param Sf
 *  Order parameters [a, b, c] for fast motion
 * @return R1
 *  Returns R1 as double
 */

double AIMF_13CR1(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * Ss, double * Sf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a double128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N, wCOCa;
	double d2x, d2y;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	double *csa;
	csa = res->csaC;

	d2x = (double) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
	d2y = (double) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
	//d2xy= (double) sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* Calculate all order parameters */
	double S2CHs, S2CSAxs, S2CSAys, S2CNs, S2CCs;
	double S2CHf, S2CSAxf, S2CSAyf, S2CNf, S2CCf;
	
	double *S2s[] = {&S2CHs, &S2CSAxs, &S2CSAys, &S2CNs, &S2CCs};
	double *S2f[] = {&S2CHf, &S2CSAxf, &S2CSAyf, &S2CNf, &S2CCf};
	
	struct Orient * As[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	
	AIMF_S2(Ss, As, S2s, 5);
	AIMF_S2(Sf, As, S2f, 5);
	
	/* N CSA relaxation contribution */
	double R1CSAx, R1CSAy, R1CSA, R1CH, R1CHr, R1CN, R1CC;
	double J1 = 0;


	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	J1 = J0_EMF(omega_13C, taus, S2CSAxs, tauf, S2CSAxf);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
	J1 = J0_EMF(omega_13C, taus, S2CSAys, tauf, S2CSAyf);
	R1CSAy = (1/15.) * d2y * J1;
	//J1 = J0_EMF_CC(omega_13C, taus, S2CSAxys, tauf, S2CSAxyf);
	//R1CSAxy = (1/15.) * d2xy * J1;
	R1CSA = R1CSAx + R1CSAy; // + 2.*R1CSAxy;
	
	

	/* N Dipolar Interactions Contributions */
	R1CH = AIMF_Dipolar_R1(omega_13C, omega_1H, taus, S2CHs, tauf, S2CHf, D_CH);
	R1CHr = AIMF_Dipolar_R1(omega_13C, omega_1H, taus, S2CNs, tauf, S2CNf, D_CHr);
	R1CN = AIMF_Dipolar_R1(omega_13C, omega_15N, taus, S2CNs, tauf, S2CNf, D_CN);
	R1CC = AIMF_Dipolar_R1(omega_13C, omega_13C - wCOCa, taus, S2CCs, tauf, S2CCf, D_CC);

	return (double) (R1CSA + R1CH + R1CHr + R1CN + R1CC)*T_DOWN;
}

/**
 * Calculates R2 relaxation rate for 13C nucleus under anisotropic model free.

 * @param res
 *  Pointer to residue being considered
 * @param relax
 *  Pointer to relaxation measurement being modelled
 * @param taus
 *  Correlation time of slow motion
 * @param tauf
 *  Correlation time of fast motion
 * @param Ss
 *  Order parameters [a, b, c] for slow motion
 * @param sigf
 *  Order parameters [a, b, c] for fast motion
 * @return R2
 *  Returns R2 as double
 */

double AIMF_13CR2(struct Residue *res, struct Relaxation* relax, double taus, double tauf, double * Ss, double * Sf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz
	double w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a double128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N, wCOCa;
	double d2x, d2y;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	double *csa;
	csa = res->csaC;

	d2x = (double) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
	d2y = (double) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
	//d2xy= (double) sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	double S2CHs, S2CSAxs, S2CSAys, S2CNs, S2CCs;
	double S2CHf, S2CSAxf, S2CSAyf, S2CNf, S2CCf;

	double *S2s[] = {&S2CHs, &S2CSAxs, &S2CSAys, &S2CNs, &S2CCs};
	double *S2f[] = {&S2CHf, &S2CSAxf, &S2CSAyf, &S2CNf, &S2CCf};
	
	struct Orient * As[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	
	AIMF_S2(Ss, As, S2s, 5);
	AIMF_S2(Sf, As, S2f, 5);
	

	/* CSA relaxation contribution */
	double R2CSAx, R2CSAy, R2CSA, R2CH, R2CHr, R2CN, R2CC;
	
	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	R2CSAx = AIMF_CSA_R2(omega_13C, w1, wr, taus, S2CSAxs, tauf, S2CSAxf, d2x, J0_EMF);
	R2CSAy = AIMF_CSA_R2(omega_13C, w1, wr, taus, S2CSAys, tauf, S2CSAyf, d2y, J0_EMF);
	//R2CSAxy = AIMF_CSA_R2(omega_13C, w1, wr, taus, S2CSAxys, tauf, S2CSAxyf, d2xy, J0_EMF_CC);
	R2CSA = R2CSAx + R2CSAy;// + 2.*R2CSAxy;

	/* N Dipolar Interactions Contributions */
	R2CH = AIMF_Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, S2CHs, tauf, S2CHf, D_CH);
	R2CHr = AIMF_Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, S2CNs, tauf, S2CNf, D_CHr);
	R2CN = AIMF_Dipolar_R2(omega_13C, omega_15N, w1, wr, taus, S2CNs, tauf, S2CNf, D_CN);
	R2CC = AIMF_Dipolar_R2(omega_13C, omega_13C - wCOCa, w1, wr, taus, S2CCs, tauf, S2CCf, D_CC);
	return (double) ((R2CSA + R2CH + R2CHr + R2CN + R2CC)*(double)T_DOWN);
}


