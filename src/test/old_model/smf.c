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
 * smf.c
 *
 * Implements simple model free.
 */

#include <stdio.h>
#include "datatypes.h"
#include "smf.h"

/**
 * Calculates spectral density function for given frequency (omega) according to simple model free analysis.
 * Implements J(w) = (1 - S2) tau / (1 + (w tau)^2) as from Lipari1982 eq 35 (assuming no overall tumbling)
 * and the factor of 2/5 has been taken into the summation (see below).
 */
Decimal J0_SMF(Decimal omega, Decimal tau, Decimal S2) {
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
Decimal SMF_Dipolar_R1(Decimal omega_X, Decimal omega_Y, Decimal d, Decimal tau, Decimal S2) {
	// omega_X is the measured nuclei (eg 15N, 13C), Y is the other (eg 1H)
	Decimal Jcomp = 0;
	Jcomp += J0_SMF(omega_Y - omega_X, tau, S2);
	Jcomp += 3. * J0_SMF(omega_X, tau, S2);
	Jcomp += 6. * J0_SMF(omega_Y + omega_X, tau, S2);
	return (Decimal) (0.1 * d * d * Jcomp);
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
Decimal SMF_Dipolar_R2(Decimal omega_X, Decimal omega_Y, Decimal d, Decimal tau, Decimal S2, Decimal J0sum) {
	// all args should've been T_DOWN'ed
	Decimal JD = J0sum + 3 * J0_SMF(omega_X, tau, S2);
	JD += J0_SMF(omega_Y - omega_X, tau, S2);
	JD += 6 * J0_SMF(omega_Y, tau, S2);
	JD += 6 * J0_SMF(omega_Y + omega_X, tau, S2);
	return (Decimal) ((1/20.) * d * d * JD);
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
 *  Returns R1 as Decimal
 */
Decimal SMF_R1(struct Residue *res, struct Relaxation* relax, Decimal tau, Decimal S2, unsigned int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */

	Decimal field = relax->field * 1000000; // conversion to Hz
	//tau = tau * T_UP;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Constants are in datatypes.c */

	Decimal omega_1H = T_DOWN * 2 * M_PI * field;
	Decimal omega_15N = T_DOWN * 2 * M_PI * field / 9.869683408806043;
	Decimal omega_13C = T_DOWN * 2 * M_PI * field / 3.976489314034722;
	Decimal omega_L;
	Decimal d2tot;
	Decimal wCOCa = 120 * omega_13C * 0.000001;
	Decimal *csa;
	Decimal R1D = 0, R1CSA;
	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		//d = -D_NH;
		// Decimal SMF_Dipolar_R1(Decimal omega_X, Decimal omega_Y, Decimal d, Decimal tau, Decimal S2) {
		R1D += SMF_Dipolar_R1(omega_15N, omega_1H, -D_NH, tau, S2); // N-H
		R1D += SMF_Dipolar_R1(omega_15N, omega_1H, -D_NHr, tau, S2); // N-Hr
		R1D += SMF_Dipolar_R1(omega_15N, omega_13C, -D_CN, tau, S2); // N-C
		R1D += SMF_Dipolar_R1(omega_15N, omega_13C, -D_NCA, tau, S2); // N-Ca
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		R1D += SMF_Dipolar_R1(omega_13C, omega_1H, -D_CH, tau, S2); // C-H
		R1D += SMF_Dipolar_R1(omega_13C, omega_1H, -D_CHr, tau, S2); // C-Hr
		R1D += SMF_Dipolar_R1(omega_13C, omega_15N, -D_CN, tau, S2); // C-N
		R1D += SMF_Dipolar_R1(omega_13C, omega_13C - wCOCa, -D_CCAc, tau, S2); // C-C
        R1D += SMF_Dipolar_R1(omega_13C, omega_13C - wCOCa, -D_CCAp, tau, S2); // C-C
        //d = -D_CH;
	} else {
		ERROR("Incorrect mode.");
		return -1;
	}

	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	omega_L *= T_DOWN;
	
	R1CSA = (2/15.) * d2tot * J0_SMF(omega_L, tau, S2);
	
	R1D = R1D * T_DOWN;
	R1CSA = R1CSA * T_DOWN;
	return (Decimal) R1D + (Decimal) R1CSA;
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
 *  Returns R2 as Decimal
 */
Decimal SMF_R2(struct Residue *res, struct Relaxation* relax, Decimal tau, Decimal S2, unsigned int mode) {

	Decimal field = relax->field * 1000000; // conversion to Hz

	//tau *= T_UP;
	
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. */

	Decimal omega_1H = T_DOWN * 2 * M_PI * field;
	Decimal omega_15N = T_DOWN * 2 * M_PI * field / 9.869683408806043;
	Decimal omega_13C = T_DOWN * 2 * M_PI * field / 3.976489314034722;
	Decimal omega_L;
	Decimal d2tot;
	Decimal *csa;
	Decimal R2D = 0, R2CSA;
	Decimal wCOCa = 120 * omega_13C * 0.000001;
	Decimal w1 = relax->w1;
	Decimal wr = relax->wr;
	w1 *= T_DOWN;
	wr *= T_DOWN;
	
	Decimal J0sum = 0;
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 - 2 * wr), tau, S2);
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 + 2 * wr), tau, S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 - wr), tau, S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 + wr), tau, S2);
	
	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		//d = -D_NH;
		R2D += SMF_Dipolar_R2(omega_15N, omega_1H, -D_NH, tau, S2, J0sum); // N-H
		R2D += SMF_Dipolar_R2(omega_15N, omega_1H, -D_NHr, tau, S2, J0sum); // N-Hr
		R2D += SMF_Dipolar_R2(omega_15N, omega_13C, -D_CN, tau, S2, J0sum); // N-C
		R2D += SMF_Dipolar_R2(omega_15N, omega_13C, -D_NCA, tau, S2, J0sum); // N-Ca
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		R2D += SMF_Dipolar_R2(omega_13C, omega_1H, -D_CH, tau, S2, J0sum); // C-H
		R2D += SMF_Dipolar_R2(omega_13C, omega_1H, -D_CHr, tau, S2, J0sum); // C-Hr
		R2D += SMF_Dipolar_R2(omega_13C, omega_15N, -D_CN, tau, S2, J0sum); // C-N
		R2D += SMF_Dipolar_R2(omega_13C, omega_13C - wCOCa, -D_CCAp, tau, S2, J0sum); // C-C
        R2D += SMF_Dipolar_R2(omega_13C, omega_13C - wCOCa, -D_CCAc, tau, S2, J0sum); // C-C

		//d = -D_CH;
	} else {
		ERROR("Incorrect mode.");
		return -1;
	}
	

	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);
	
	omega_L *= T_DOWN;

	R2CSA = (1/45.) * d2tot * (J0sum + 3 * J0_SMF(omega_L, tau, S2));
	R2D *= T_DOWN;
	R2CSA *= T_DOWN;
	return (Decimal) R2D + (Decimal) R2CSA;
}
