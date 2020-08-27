/**
 * @file models.c
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
  

long double sq(long double x) {
	return x * x;
}

int sq_i(int x) {
	return x * x;
}


  
/** 
 * Calculates R1 contribution under GAF model for Dipolar interaction between two atoms.
 * omega_obs is the frequency of the observed nucleus, omega_neigh is the frequency of the other atoms.\n
 * D is the dipolar coupling in rad/s.
 * Equation as 3.28, McConnell1986, only for solid state spins.
 * Kurbanov2011, eq 9 (factor of 1/10 arising from (2/5.) in spectral density and 1/4. in dipolar.
 */
inline double GAF_Dipolar_R1(long double omega_obs, long double omega_neigh, long double taus, long double S2s, long double tauf, long double S2f, long double D) {
	return (0.1) * sq(D) * ( \
		(J0_EMF(omega_neigh - omega_obs, taus, S2s, tauf, S2f)) \
		 + 3 * (J0_EMF(omega_obs, taus, S2s, tauf, S2f)) \
		  + 6 * (J0_EMF(omega_neigh + omega_obs, taus, S2s, tauf, S2f)) \
	  );
}
    

/**
 * Calculates the R2 contribution of CSA of nuclei with frequency omega.
 * D2 is the squared D22x/D22y/D22xy variable (dependent on CSA)
 * Kurbanov2011, eqs 18-20.
 */
inline double GAF_CSA_R2(long double omega, float w1, float wr, long double taus, long double S2s, long double tauf, long double S2f, double D2) {
	return ( \
  		(1/45.) * (double) D2 * ( \
  			(2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f) +\
  			(2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f) + \
  			(3.)   * J0_EMF(omega, taus, S2s, tauf, S2f)\
  		)\
  	);
}

/** 
 * Calculates R2 contribution of dipolar interaction between two atoms as in GAF_Dipolar_R1
 * Kurbanov2011, for the case where theta_p = 90 such that sin(theta_p) = 1 and so
 * R1p = 0.5*R1 + R1delta
 */
inline double GAF_Dipolar_R2(long double omega_obs, long double omega_neigh, float w1, float wr, long double taus, long double S2s, long double tauf, long double S2f, double D) {
	return (\
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
double SMF_R1(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */

	long double field = relax->field * 1000000; // conversion to Hz
	//tau = tau * T_UP;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Constants are in datatypes.c */

	long double omega_1H = 2 * M_PI * field;
	long double omega_L;
	long double d = 0, d2tot;
	float *csa;

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
double SMF_R2(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, int mode) {

	long double field = relax->field * 1000000; // conversion to Hz

	//tau *= T_UP;
	
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. */

	long double omega_1H = 2 * M_PI * field;

	long double omega_L;
	long double d, d2tot;
	float *csa;

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
	float *csa;

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
	float *csa;

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

/**
 * Calculates the order parameter relating to two orientations A and B undergoing anisotropic axial fluctuations of magnitude sig.
 * This function parallelises the calculation such that the loop over l, m, mp, k, kp only needs to occur once per motion - also
 * meaning that the calculation of the exponential only needs to be done once. Checks for conjugation are made to reduce the need
 * to calculate complex powers only to the cases where this is necessary.\n
 * Lienin1998, eq 6
 * @param sig
 *  Array containing [alpha, beta, gamma] deflection angles
 * @param A
 *  Array containing pointers to orientation vectors for 'A' for each S2 being calculated.
 * @param B
 *  Array containing pointers to orientation vectors for 'B' for each S2 being calculated.
 * @param S2
 *  Array containing pointers to S2 values - this is where output will be put.
 * @param length
 *  Number of S2 values being calculated (eg length of S2)
 * @param mode
 *  MODE_REAL or MODE_IMAG depending on which result is being returned.
 * @return Returns 1 in all cases.
 */
int GAF_S2(long double sig[3], struct Orient ** A, struct Orient ** B, double * S2[], int length,  int mode) {
	int l, m, mp, k, kp, i;
	double complex * Amp = (double complex *) malloc(sizeof(double complex) * length);
	double complex temp, ttemp;
	for (i = 0; i < length; i++) {
		Amp[i] = 0;
	}
	
	int ksqsum;
	int lsqsum;
	int msqsum;
	double lexp, kexp, mexp;
	for (l = -2; l <= 2; l++) {
		lsqsum = sq_i(l);
		lexp = -(sq(sig[1]) * lsqsum);
		for (m = -2; m <= 2; m++) {
			for (mp = -2; mp <= 2; mp++) {
				msqsum = sq_i(m) + sq_i(mp);
				mexp = -(sq(sig[2]) * msqsum/2.);
				for (k = -2; k <= 2; k++) {
					for (kp = -2; kp <= 2; kp++) {
						ksqsum = sq_i(k) + sq_i(kp);
						kexp = -(sq(sig[0]) * ksqsum / 2.);
						//if (A->Y2[m+2] == A->Y2c[m+2] && B->Y2[mp+2] == B->Y2c[mp+2]) {
						temp = 1;
						
						temp *= exp(lexp + kexp + mexp);
						temp *= Dwig[k+2][l+2] * Dwig[kp+2][l+2] * Dwig[m+2][k+2] * Dwig[mp+2][kp+2];
						for (i = 0; i < length; i++) {
							ttemp = temp;
							if (cimag(A[i]->Y2[m+2]) < 0.00001 && cimag(B[i]->Y2[mp+2]) < 0.00001 && mode==MODE_REAL) {
								/* In this case then Y2 and Y2c are real numbers and therefore the last step of this is real.
								* If mode is MODE_REAL, then, we may safely ignore any case in which
								*  cpowl(-1*I, k-kp)
								* is imaginary (as this will give only an imaginary component to the Amp).
								* This is the case for any (k-kp)%2 != 0.
								* In this case, if (k-kp)%4 == 0 the cpowl function gives 1, else -1.
								*/
								if ((k-kp)%2 != 0)
									continue;
								if ((k - kp) % 4 == 0)
									ttemp *= 1;
								else
									ttemp *= -1.;
							} else {
								/* Other wise we have to go the long route... */
								ttemp *= cpow(-1 * I, k - kp);
							}
							
							
							
							ttemp *= A[i]->Y2[m+2] * B[i]->Y2c[mp+2];
							Amp[i] += ttemp;
						}
						//temp *= A->Y2[m+2] * B->Y2c[mp+2];
						//Amp += temp;
					}
				}
			}
		}
	}
	//Amp *= (4 * M_PI / 5.);
	for (i = 0;i < length; i++) {
		Amp[i] = Amp[i] * (4 * M_PI / 5.);
		switch (mode) {
			case MODE_REAL: *S2[i] = creal(Amp[i]); break;
			case MODE_IMAG: *S2[i] = creal(Amp[i]); break;
			default: break;
		}
	}
	free(Amp);
	return 1;
}





/* This has not been directly tested against the MATLAB model as I don't really
 * have an equivalent to it set up (NHCO_3GAFSquaredEa is likely the closest
 * but I don't have the setup around it for the Legendre polynomials) */
 /* WARNING: NCSAxy IS OMITTED */
 
/**
 * Calculates R1 relaxation rate for 15N nucleus under gaussian axial fluctuations.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @warning I'm unsure of the validity of NCSAxy.
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
double GAF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
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

	float *csa;
	csa = res->csaN;

	d2x = sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* Calculate all order parameters */
	double S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	double S2NHf, S2NCSAxf, S2NCSAyf, S2NCSAxyf, S2CNf, S2CaNf;
	
	double *S2s[] = {&S2NHs, &S2NCSAxs, &S2NCSAys, &S2NCSAxys, &S2CNs, &S2CaNs};
	double *S2f[] = {&S2NHf, &S2NCSAxf, &S2NCSAyf, &S2NCSAxyf, &S2CNf, &S2CaNf};
	
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	struct Orient * Bs[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	GAF_S2(sigf, As, Bs, S2f, 6, MODE_REAL);
	

	/* N CSA relaxation contribution */
	double R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1NH, R1NHr, R1CN, R1CaN;
	long double J1 = 0;
	

	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	

	J1 = J0_EMF(omega_15N, taus, S2NCSAxs, tauf, S2NCSAxf);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
	J1 = J0_EMF(omega_15N, taus, S2NCSAys, tauf, S2NCSAyf);
	R1CSAy = (1/15.) * d2y * J1;
	J1 = J0_EMF(omega_15N, taus, S2NCSAxys, tauf, S2NCSAxyf);
	R1CSAxy = (2/15.) * d2xy * J1;
	
	/** Eq A30, Bremi1997 */
	R1CSA = R1CSAx + R1CSAy + R1CSAxy;

	printf("ts: %Le, S2xys: %f, tf: %Le, S2xyf: %f\n\tR1CSAxy = %f\n", taus*T_DOWN, S2NCSAxys, tauf*T_DOWN, S2NCSAxyf, R1CSAxy*T_DOWN);
	printf("\n\tR1CSA (%f) = R1CSAx (%f) + R1CSAy (%f) - R1CSAxy (%f)\n", R1CSA*T_DOWN, R1CSAx*T_DOWN, R1CSAy*T_DOWN, R1CSAxy*T_DOWN);



	/* N Dipolar Interactions Contributions */
	R1NH = GAF_Dipolar_R1(omega_15N, omega_1H, taus, S2NHs, tauf, S2NHf, D_NH);
	R1NHr = GAF_Dipolar_R1(omega_15N, omega_1H, taus, S2CNs, tauf, S2CNf, D_HNr);
	R1CN = GAF_Dipolar_R1(omega_15N, omega_13C, taus, S2CNs, tauf, S2CNf, D_CN);
	R1CaN = GAF_Dipolar_R1(omega_15N, omega_13C, taus, S2CaNs, tauf, S2CaNf, D_CaN);

	return (R1CSA + R1NH + R1NHr + R1CN + R1CaN) * T_DOWN;
}

/**
 * Calculates R2 relaxation rate for 15N nucleus under gaussian axial fluctuations.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @warning I'm unsure of the validity of NCSAxy.
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
 * @return R2
 *  Returns R2 as double
 */
double GAF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz
	float w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N;
	double d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;

	float *csa;
	csa = res->csaN;

	d2x = sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	double S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	double S2NHf, S2NCSAxf, S2NCSAyf, S2NCSAxyf, S2CNf, S2CaNf;
	
	double *S2s[] = {&S2NHs, &S2NCSAxs, &S2NCSAys, &S2NCSAxys, &S2CNs, &S2CaNs};
	double *S2f[] = {&S2NHf, &S2NCSAxf, &S2NCSAyf, &S2NCSAxyf, &S2CNf, &S2CaNf};
	
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	struct Orient * Bs[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	GAF_S2(sigs, As, Bs, S2f, 6, MODE_REAL);

	/* N CSA relaxation contribution */
	double R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2NH, R2NHr, R2CN, R2CaN;

	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	R2CSAx = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAxs, tauf, S2NCSAxf, d2x);
	R2CSAy = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAys, tauf, S2NCSAyf, d2y);
	R2CSAxy = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAxys, tauf, S2NCSAxyf, d2xy);
	/**
	 * I'm unsure if this should be R2CSAx + R2CSAy + R2CSAxy or - R2CSAxy. 
	 * In the original MATLAB scripts, it is - R2CSAxy; in NHCO_3GAFSquaredEa, the form of R2CSAxy
	 * is really odd. In NH3DEglobalGAF it is also - R2CSAxy but the form is the same.
	 * In Lienin1998's SI, eq 5 and 6 it is +R2CSAxy - same for Bremi1997 eq A30
	 */
	R2CSA = R2CSAx + R2CSAy + R2CSAxy;

	/* N Dipolar Interactions Contributions */
	R2NH = GAF_Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, S2NHs, tauf, S2NHf, D_NH);
	R2NHr = GAF_Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, S2CNs, tauf, S2CNf, D_HNr);
	R2CN = GAF_Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, S2CNs, tauf, S2CNf, D_CN);
	R2CaN = GAF_Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, S2CaNs, tauf, S2CaNf, D_CaN);
	return (R2CSA + R2NH + R2NHr + R2CN + R2CaN)*T_DOWN;
}



/* This has not been directly tested against the MATLAB model as I don't really
 * have an equivalent to it set up (NHCO_3GAFSquaredEa is likely the closest
 * but I don't have the setup around it for the Legendre polynomials) */
 /* WARNING: CSAxy IS OMITTED */
 
 /**
 * Calculates R1 relaxation rate for 13C nucleus under gaussian axial fluctuations.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @warning Has not been tested on 13C data
 * @warning I'm unsure of the validity of CSAxy.
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

double GAF_13CR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
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

	float *csa;
	csa = res->csaC;

	d2x = sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
	d2y = sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
	d2xy= sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* Calculate all order parameters */
	double S2CHs, S2CSAxs, S2CSAys, S2CSAxys, S2CNs, S2CCs;
	double S2CHf, S2CSAxf, S2CSAyf, S2CSAxyf, S2CNf, S2CCf;
	
	double *S2s[] = {&S2CHs, &S2CSAxs, &S2CSAys, &S2CSAxys, &S2CNs, &S2CCs};
	double *S2f[] = {&S2CHf, &S2CSAxf, &S2CSAyf, &S2CSAxyf, &S2CNf, &S2CCf};
	
	struct Orient * As[] = {&(res->orients[OR_CH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_CCAc])};
	struct Orient * Bs[] = {&(res->orients[OR_CH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_CCAc])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	GAF_S2(sigs, As, Bs, S2f, 6, MODE_REAL);
	
	/* N CSA relaxation contribution */
	double R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1CH, R1CHr, R1CN, R1CC;
	long double J1 = 0;


	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	J1 = J0_EMF(omega_13C, taus, S2CSAxs, tauf, S2CSAxf);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
	J1 = J0_EMF(omega_13C, taus, S2CSAys, tauf, S2CSAyf);
	R1CSAy = (1/15.) * d2y * J1;
	J1 = J0_EMF(omega_13C, taus, S2CSAxys, tauf, S2CSAxyf);
	R1CSAxy = (2/15.) * d2xy * J1;
	R1CSA = R1CSAx + R1CSAy + R1CSAxy;
	
	

	/* N Dipolar Interactions Contributions */
	R1CH = GAF_Dipolar_R1(omega_13C, omega_1H, taus, S2CHs, tauf, S2CHf, D_CH);
	R1CHr = GAF_Dipolar_R1(omega_13C, omega_1H, taus, S2CNs, tauf, S2CNf, D_CHr);
	R1CN = GAF_Dipolar_R1(omega_13C, omega_15N, taus, S2CNs, tauf, S2CNf, D_CN);
	R1CC = GAF_Dipolar_R1(omega_13C, omega_13C + wCOCa, taus, S2CCs, tauf, S2CCf, D_CC);

	return (R1CSA + R1CH + R1CHr + R1CN + R1CC)*T_DOWN;
}

/**
 * Calculates R2 relaxation rate for 13C nucleus under gaussian axial fluctuations.
 * @warning Has not been verified against MATLAB model as there isn't really an equivalent to it that I have set up.
 * @warning Has not been tested on 13C data
 * @warning I'm unsure of the validity of CSAxy.
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
 * @return R2
 *  Returns R2 as double
 */

double GAF_13CR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz
	float w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N, wCOCa;
	double d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	float *csa;
	csa = res->csaC;

	d2x = sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
	d2y = sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
	d2xy= sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	double S2CHs, S2CSAxs, S2CSAys, S2CSAxys, S2CNs, S2CCs;
	double S2CHf, S2CSAxf, S2CSAyf, S2CSAxyf, S2CNf, S2CCf;

	double *S2s[] = {&S2CHs, &S2CSAxs, &S2CSAys, &S2CSAxys, &S2CNs, &S2CCs};
	double *S2f[] = {&S2CHf, &S2CSAxf, &S2CSAyf, &S2CSAxyf, &S2CNf, &S2CCf};
	
	struct Orient * As[] = {&(res->orients[OR_CH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_CCAc])};
	struct Orient * Bs[] = {&(res->orients[OR_CH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_CCAc])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	GAF_S2(sigs, As, Bs, S2f, 6, MODE_REAL);
	

	/* CSA relaxation contribution */
	double R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2CH, R2CHr, R2CN, R2CC;
	
	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	R2CSAx = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAxs, tauf, S2CSAxf, d2x);
	R2CSAy = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAys, tauf, S2CSAyf, d2y);
	R2CSAxy = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAxys, tauf, S2CSAxyf, d2xy);
	R2CSA = R2CSAx + R2CSAy + R2CSAxy;

	/* N Dipolar Interactions Contributions */
	R2CH = GAF_Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, S2CHs, tauf, S2CHf, D_CH);
	R2CHr = GAF_Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, S2CNs, tauf, S2CNf, D_CHr);
	R2CN = GAF_Dipolar_R2(omega_13C, omega_15N, w1, wr, taus, S2CNs, tauf, S2CNf, D_CN);
	R2CC = GAF_Dipolar_R2(omega_13C, omega_13C + wCOCa, w1, wr, taus, S2CCs, tauf, S2CCf, D_CC);
	return (R2CSA + R2CH + R2CHr + R2CN + R2CC)*T_DOWN;
}
