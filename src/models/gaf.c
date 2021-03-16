/*
 * Copyright (c) 2021 Ben Tatman, University of Warwick
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the Sosftware
 * is furnished to do so, subject to the following conditions:
 *
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software
 * - Any academic work deriving from the Software shall cite [CITATION].
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUR OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * gaf.c
 *
 * Implements Gaussian Axial Fluctuation model
 */

#include <stdlib.h>
#include <complex.h>

#include "../datatypes.h"
#include "emf.h"
#include "gaf.h"
  
/** 
 * Calculates R1 contribution under GAF model for Dipolar interaction between two atoms.
 * omega_obs is the frequency of the observed nucleus, omega_neigh is the frequency of the other atoms.\n
 * D is the dipolar coupling in rad/s.
 * Equation as 3.28, McConnell1986, only for solid state spins.
 * Kurbanov2011, eq 9 (factor of 1/10 arising from (2/5.) in spectral density and 1/4. in dipolar.
 */
Decimal GAF_Dipolar_R1(Decimal omega_obs, Decimal omega_neigh, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal D) {
	Decimal q = (0.1) * sq(D) * ( \
		(J0_EMF(omega_neigh - omega_obs, taus, S2s, tauf, S2f)) \
		 + 3 * (J0_EMF(omega_obs, taus, S2s, tauf, S2f)) \
		  + 6 * (J0_EMF(omega_neigh + omega_obs, taus, S2s, tauf, S2f)) \
	  );
	return (Decimal) q;
}
    

/**
 * Calculates the R2 contribution of CSA of nuclei with frequency omega.
 * D2 is the squared D22x/D22y/D22xy variable (dependent on CSA)
 * Kurbanov2011, eqs 18-20.
 */
 //Decimal J0_EMF_CC(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f) {
 
Decimal GAF_CSA_R2(Decimal omega, \
						 Decimal w1, \
						 Decimal wr, \
						 Decimal taus, \
						 Decimal S2s, \
						 Decimal tauf, \
						 Decimal S2f, \
						 Decimal D2, \
						 Decimal  (*J_SD)(\
						 	Decimal,\
						 	Decimal, \
						 	Decimal, \
						 	Decimal, \
						 	Decimal)\
						) {
	return (Decimal) ( \
  		(1/45.) * (Decimal) D2 * ( \
  			(2/3.) * J_SD(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f) +\
  			(2/3.) * J_SD(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J_SD(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J_SD(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f) + \
  			(3.)   * J_SD(omega, taus, S2s, tauf, S2f)\
  		)\
  	);
}

/** 
 * Calculates R2 contribution of dipolar interaction between two atoms as in GAF_Dipolar_R1
 * Kurbanov2011, for the case where theta_p = 90 such that sin(theta_p) = 1 and so
 * R1p = 0.5*R1 + R1delta
 */
Decimal GAF_Dipolar_R2(Decimal omega_obs, Decimal omega_neigh, Decimal w1, Decimal wr, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal D) {
	return (Decimal) (\
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
int GAF_S2(Decimal sig[3], struct Orient ** A, struct Orient ** B, Decimal * S2[], int length,  int mode) {
	int l, m, mp, k, kp, i;
	Complex * Amp = (Complex *) malloc(sizeof(Complex) * (long unsigned int) length);
	Complex temp, ttemp;


	for (i = 0; i < length; i++) {
		Amp[i] = 0;
	}
	
	int ksqsum;
	int lsqsum;
	int msqsum;

    Decimal sqsigs[3] = {sq(sig[0]), sq(sig[1]), sq(sig[2])};

    /* sq_s is an array where sq_s[i+2] = i^2 */
    int sq_s[5] = {4, 1, 0, 1, 4};

    /* sq_d is an array where sq_d[i+2][j+2] = (i^2 + j^2) / 2
     * This speeds precalculates loop values */
    Decimal sq_d[5][5] = {
            {4.0, 2.5, 2.0, 2.5, 4.0},
            {2.5, 1.0, 0.5, 1.0, 2.5},
            {2.0, 0.5, 0.0, 0.5, 2.0},
            {2.5, 1.0, 0.5, 1.0, 2.5},
            {4.0, 2.5, 2.0, 2.5, 4.0}
    };

	Decimal lexp, kexp, mexp;

	/*
	 * The GAF order parameter equation is composed of four main parts. These are
	 *  1. Imaginary power of k-kp
	 *  2. Exponential of angular deflections with (m^2 + mp^2) etc
	 *  3. Product of Wigner matrices
	 *  4. Product of Spherical Harmonics
	 * As this is a nested sum it can be quite slow. In order to speed it up, a number of optimizations have been made.
	 *
	 *  1. If we are considering only real order parameters, then given spherical harmonics with
	 *     no imaginary component it can be shown that any imaginary part of the i^(k-kp) is ignored.
	 *     As (k-kp) is an integer, i^(k-kp) goes (0) 1, (1) i, (2) -1, (3) -i, (4) 1 and so on.
	 *     So creal(i^(k)) is 0 for any odd k, and so these cases may be ignored.
	 *     If k is divisible by 4, then this factor is 1, and if not, it is -1.
	 *     This avoids the need for calling cpow() on each iteration, which is a hard operation.
	 *     In the case that the spherical harmonics have an imaginary part, however, this must be done.
	 *  2. The component for each of alpha, beta, and gamma is done at their section in the loop. This means
	 *     that the beta component is only calculated 5 times, the gamma component 125 times, and the alpha
	 *     3125 times.
	 *     The sigma values are presquared in sqsigs. While squaring is easy, doing this so many times can slow
	 *     down the iteration. Additionally, the summed squares of the indices is precalculated (sq_s, sq_d).
	 *     It has been considered to precalculate the sigma exponentials and then take these to the power of the
	 *     sum of squares. This doesn't lead to a significant speed bump, but does make the code more complicated.
	 *  3. The product of Wigner matrices over k, kp, l, m, mp and exponentiation is done once and this is then
	 *     used for all S2 parameters being calculated.
	 *  4. No optimization is done to the products of spherical harmonics.
	 */

	for (l = -2; l <= 2; l++) {
		//lsqsum = sq_i(l);
		//lexp = (Decimal) -(sq(sig[1]) * lsqsum);
		lexp = -(sqsigs[1] * sq_s[l+2]);
		for (m = -2; m <= 2; m++) {
			for (mp = -2; mp <= 2; mp++) {
				//msqsum = sq_i(m) + sq_i(mp);
				//mexp = (Decimal) -(sq(sig[2]) * msqsum/2.);
				mexp = -(sqsigs[2] * sq_d[m+2][mp+2]);
				for (k = -2; k <= 2; k++) {
					for (kp = -2; kp <= 2; kp++) {
						//ksqsum = sq_i(k) + sq_i(kp);
						//kexp = (Decimal) -(sq(sig[0]) * ksqsum / 2.);
						kexp = -(sqsigs[0] * sq_d[k+2][kp+2]);
						//if (A->Y2[m+2] == A->Y2c[m+2] && B->Y2[mp+2] == B->Y2c[mp+2]) {
						temp = 1;
						
						temp *= exp(lexp + kexp + mexp);
						temp *= Dwig[k+2][l+2] * Dwig[kp+2][l+2] * Dwig[m+2][k+2] * Dwig[mp+2][kp+2];
						for (i = 0; i < length; i++) {
							ttemp = temp;
							if (cimag(A[i]->Y2[m+2]) < 0.00001 && cimag(B[i]->Y2[mp+2]) < 0.00001 && mode==MODE_REAL) {
								/* In this case then Y2 and Y2c are real numbers and therefore the last step of this is real.
								* If mode is MODE_REAL, then, we may safely ignore any case in which
								*  cpow(-1*I, k-kp)
								* is imaginary (as this will give only an imaginary component to the Amp).
								* This is the case for any (k-kp)%2 != 0.
								* In this case, if (k-kp)%4 == 0 the cpow function gives 1, else -1.
								*/
								if ((k-kp)%2 != 0)
									continue;
								if ((k - kp) % 4 == 0)
									ttemp *= 1;
								else
									ttemp *= -1.;
							} else {
								//printf("Going long route...\n");
								/* Other wise we have to go the long route... */
								ttemp *= cpow(-1 * I, k - kp);
							}
							
							
							ttemp *= A[i]->Y2[m+2] * B[i]->Y2c[mp+2];

                            Amp[i] += ttemp;
							//printf("Amp %d %lf\n", i, Amp[i]);
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
			case MODE_REAL: *S2[i] = creal((Complex) Amp[i]); break;
			case MODE_IMAG: *S2[i] = cimag((Complex) Amp[i]); break;
			default: break;
		}
	}
	free(Amp);
	return 1;
}





/**
 * Calculates R1 relaxation rate for 15N nucleus under gaussian axial fluctuations.
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
 * @param sigf
 *  Deflection angles [alpha, beta, gamma] for fast motion
 * @return R1
 *  Returns R1 as Decimal
 */
Decimal GAF_15NR1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	Decimal field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	Decimal omega_1H = 2 * M_PI * field;
	Decimal omega_13C, omega_15N;
	Decimal d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;

	Decimal *csa;
	csa = res->csaN;

	d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= (Decimal) sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* Calculate all order parameters */
	Decimal S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	Decimal S2NHf, S2NCSAxf, S2NCSAyf, S2NCSAxyf, S2CNf, S2CaNf;
	
	Decimal *S2s[] = {&S2NHs, &S2NCSAxs, &S2NCSAys, &S2NCSAxys, &S2CNs, &S2CaNs};
	Decimal *S2f[] = {&S2NHf, &S2NCSAxf, &S2NCSAyf, &S2NCSAxyf, &S2CNf, &S2CaNf};
	
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	struct Orient * Bs[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	GAF_S2(sigf, As, Bs, S2f, 6, MODE_REAL);
	

	/* N CSA relaxation contribution */
	Decimal R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1NH, R1NHr, R1CN, R1CaN;
	Decimal J1;
	

	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	

	J1 = J0_EMF(omega_15N, taus, S2NCSAxs, tauf, S2NCSAxf);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
	J1 = J0_EMF(omega_15N, taus, S2NCSAys, tauf, S2NCSAyf);
	R1CSAy = (1/15.) * d2y * J1;
	J1 = J0_EMF_CC(omega_15N, taus, S2NCSAxys, tauf, S2NCSAxyf);
	R1CSAxy = (1/15.) * d2xy * J1;
	
	/** Eq A30, Bremi1997 */
	R1CSA = R1CSAx + R1CSAy + 2.*R1CSAxy;
	



	/* N Dipolar Interactions Contributions */
	R1NH = GAF_Dipolar_R1(omega_15N, omega_1H, taus, S2NHs, tauf, S2NHf, D_NH);
	R1NHr = GAF_Dipolar_R1(omega_15N, omega_1H, taus, S2CNs, tauf, S2CNf, D_HNr);
	R1CN = GAF_Dipolar_R1(omega_15N, omega_13C, taus, S2CNs, tauf, S2CNf, D_CN);
	R1CaN = GAF_Dipolar_R1(omega_15N, omega_13C, taus, S2CaNs, tauf, S2CaNf, D_CaN);

	return (Decimal) (R1CSA + R1NH + R1NHr + R1CN + R1CaN) * T_DOWN;
}

/**
 * Calculates R2 relaxation rate for 15N nucleus under gaussian axial fluctuations.
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
 * @return R2
 *  Returns R2 as Decimal
 */
Decimal GAF_15NR2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	Decimal field = relax->field * 1000000; // conversion to Hz
	Decimal w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	Decimal omega_1H = 2 * M_PI * field;
	Decimal omega_13C, omega_15N;
	Decimal d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;

	Decimal *csa;
	csa = res->csaN;

	d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= (Decimal) sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	Decimal S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	Decimal S2NHf, S2NCSAxf, S2NCSAyf, S2NCSAxyf, S2CNf, S2CaNf;
	
	Decimal *S2s[] = {&S2NHs, &S2NCSAxs, &S2NCSAys, &S2NCSAxys, &S2CNs, &S2CaNs};
	Decimal *S2f[] = {&S2NHf, &S2NCSAxf, &S2NCSAyf, &S2NCSAxyf, &S2CNf, &S2CaNf};
	
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	struct Orient * Bs[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	GAF_S2(sigf, As, Bs, S2f, 6, MODE_REAL);

	/* N CSA relaxation contribution */
	Decimal R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2NH, R2NHr, R2CN, R2CaN;

	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	R2CSAx = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAxs, tauf, S2NCSAxf, d2x, J0_EMF);
	R2CSAy = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAys, tauf, S2NCSAyf, d2y, J0_EMF);
	R2CSAxy = GAF_CSA_R2(omega_15N, w1, wr, taus, S2NCSAxys, tauf, S2NCSAxyf, d2xy, J0_EMF_CC);
	/**
	 * Equation for spectral density for xy taken from Manuscript
	 */ 
	R2CSA = R2CSAx + R2CSAy + 2.* R2CSAxy;

	/* N Dipolar Interactions Contributions */
	R2NH = GAF_Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, S2NHs, tauf, S2NHf, D_NH);
	R2NHr = GAF_Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, S2CNs, tauf, S2CNf, D_HNr);
	R2CN = GAF_Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, S2CNs, tauf, S2CNf, D_CN);
	R2CaN = GAF_Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, S2CaNs, tauf, S2CaNf, D_CaN);
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
 *  Returns R1 as Decimal
 */

Decimal GAF_13CR1(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	Decimal field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	Decimal omega_1H = 2 * M_PI * field;
	Decimal omega_13C, omega_15N, wCOCa;
	Decimal d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	Decimal *csa;
	csa = res->csaC;

	d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
	d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
	d2xy= (Decimal) sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* Calculate all order parameters */
	Decimal S2CHs, S2CSAxs, S2CSAys, S2CSAxys, S2CNs, S2CCs;
	Decimal S2CHf, S2CSAxf, S2CSAyf, S2CSAxyf, S2CNf, S2CCf;
	
	Decimal *S2s[] = {&S2CHs, &S2CSAxs, &S2CSAys, &S2CSAxys, &S2CNs, &S2CCs};
	Decimal *S2f[] = {&S2CHf, &S2CSAxf, &S2CSAyf, &S2CSAxyf, &S2CNf, &S2CCf};
	
	struct Orient * As[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	struct Orient * Bs[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	GAF_S2(sigf, As, Bs, S2f, 6, MODE_REAL);
	
	/* N CSA relaxation contribution */
	Decimal R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1CH, R1CHr, R1CN, R1CC;
	Decimal J1;


	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	J1 = J0_EMF(omega_13C, taus, S2CSAxs, tauf, S2CSAxf);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
	J1 = J0_EMF(omega_13C, taus, S2CSAys, tauf, S2CSAyf);
	R1CSAy = (1/15.) * d2y * J1;
	J1 = J0_EMF_CC(omega_13C, taus, S2CSAxys, tauf, S2CSAxyf);
	R1CSAxy = (1/15.) * d2xy * J1;
	R1CSA = R1CSAx + R1CSAy + 2.*R1CSAxy;
	
	

	/* N Dipolar Interactions Contributions */
	R1CH = GAF_Dipolar_R1(omega_13C, omega_1H, taus, S2CHs, tauf, S2CHf, D_CH);
	R1CHr = GAF_Dipolar_R1(omega_13C, omega_1H, taus, S2CNs, tauf, S2CNf, D_CHr);
	R1CN = GAF_Dipolar_R1(omega_13C, omega_15N, taus, S2CNs, tauf, S2CNf, D_CN);
	R1CC = GAF_Dipolar_R1(omega_13C, omega_13C - wCOCa, taus, S2CCs, tauf, S2CCf, D_CC);

	return (Decimal) (R1CSA + R1CH + R1CHr + R1CN + R1CC)*T_DOWN;
}

/**
 * Calculates R2 relaxation rate for 13C nucleus under gaussian axial fluctuations.
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
 * @return R2
 *  Returns R2 as Decimal
 */

Decimal GAF_13CR2(struct Residue *res, struct Relaxation* relax, Decimal taus, Decimal tauf, Decimal * sigs, Decimal * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	Decimal field = relax->field * 1000000; // conversion to Hz
	Decimal w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	Decimal omega_1H = 2 * M_PI * field;
	Decimal omega_13C, omega_15N, wCOCa;
	Decimal d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	Decimal *csa;
	csa = res->csaC;

	d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
	d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
	d2xy= (Decimal) sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	Decimal S2CHs, S2CSAxs, S2CSAys, S2CSAxys, S2CNs, S2CCs;
	Decimal S2CHf, S2CSAxf, S2CSAyf, S2CSAxyf, S2CNf, S2CCf;

	Decimal *S2s[] = {&S2CHs, &S2CSAxs, &S2CSAys, &S2CSAxys, &S2CNs, &S2CCs};
	Decimal *S2f[] = {&S2CHf, &S2CSAxf, &S2CSAyf, &S2CSAxyf, &S2CNf, &S2CCf};
	
	struct Orient * As[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	struct Orient * Bs[] = {&(res->orients[OR_CNH]), &(res->orients[OR_CCSAxx]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_CCAp])};
	
	GAF_S2(sigs, As, Bs, S2s, 6, MODE_REAL);
	GAF_S2(sigf, As, Bs, S2f, 6, MODE_REAL);
	

	/* CSA relaxation contribution */
	Decimal R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2CH, R2CHr, R2CN, R2CC;
	
	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	
	R2CSAx = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAxs, tauf, S2CSAxf, d2x, J0_EMF);
	R2CSAy = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAys, tauf, S2CSAyf, d2y, J0_EMF);
	R2CSAxy = GAF_CSA_R2(omega_13C, w1, wr, taus, S2CSAxys, tauf, S2CSAxyf, d2xy, J0_EMF_CC);
	R2CSA = R2CSAx + R2CSAy + 2.*R2CSAxy;

	/* N Dipolar Interactions Contributions */
	R2CH = GAF_Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, S2CHs, tauf, S2CHf, D_CH);
	R2CHr = GAF_Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, S2CNs, tauf, S2CNf, D_CHr);
	R2CN = GAF_Dipolar_R2(omega_13C, omega_15N, w1, wr, taus, S2CNs, tauf, S2CNf, D_CN);
	R2CC = GAF_Dipolar_R2(omega_13C, omega_13C - wCOCa, w1, wr, taus, S2CCs, tauf, S2CCf, D_CC);
	return (Decimal) ((R2CSA + R2CH + R2CHr + R2CN + R2CC)*(Decimal)T_DOWN);
}
