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
 * datatypes.c
 *
 * Defines variables and parameters used for the various models. Defines parameter limits
 * and Wigner matrices, and provides functions to generate parameters and obtain uniformly
 * distributed random numbers.
 */

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "datatypes.h"

time_t start_time;
Decimal Dwig[5][5];

Decimal sq(Decimal x) {
	return x * x;
}

Decimal temp_tau(Decimal tau0, Decimal Ea, Decimal temp) {
	Decimal tfactor = (300. - temp) / (temp * 300.);
	return tau0 * exp((Ea / RYD) * tfactor); // tau0 now at 300 K.
}

/**
 * Populates the second order spherical harmonics for orientation vector passed as pointer.
 * @param or
 *  Orientation vector
 * @return NULL
 */
void calculate_Y2(struct Orient * or) {
	/* or->Y2[0] is Y2 -2
	   or->Y2[1] is Y2 -1
	   or->Y2[2] is Y2 0
	   or->Y2[3] is Y2 1
	   or->Y2[4] is Y2 2 */

	/* as given https://mathworld.wolfram.com/SphericalHarmonic.html*/
	/* Y2[0] is Y2 -2*/

	or->Y2[4] = (1/4.) * (sqrt(15. / (2 * M_PI))) * (pow(sin(or->theta), 2.)) * cexp(2 * or->phi * I);
	or->Y2[3] = (-1/2.) * (sqrt(15. / (2 * M_PI))) * sin(or->theta) * cos(or->theta) * cexp(or->phi * I);
	or->Y2[2] = (1/4.) * (sqrt(5. / M_PI)) * (3 * pow(cos(or->theta), 2) - 1);
	or->Y2[1] = (1/2.) * (sqrt(15. / (2 * M_PI))) * sin(or->theta) * cos(or->theta) * cexp(-or->phi * I);
	or->Y2[0] = (1/4.) * (sqrt(15. / (2 * M_PI))) * (pow(sin(or->theta), 2.)) * cexp(-2 * or->phi * I);

	or->Y2c[0] = conj(or->Y2[0]);
	or->Y2c[1] = conj(or->Y2[1]);
	or->Y2c[2] = conj(or->Y2[2]);
	or->Y2c[3] = conj(or->Y2[3]);
	or->Y2c[4] = conj(or->Y2[4]);

}

/**
 * Initialises Wigner D matrix in global space. 
 */
void initialise_dwig(Decimal angle, Decimal Dw[5][5]) {
	// beta taken as pi/2.
	/* It could be taken as an argument but hopefully
	 * doing it this way will mean the compiler will fill it out
	 * before execution, saving a little time. */
	// 0 -2
	// 1 -1
	// 2  0
	// 3  1
	// 4  2

	// d[m' + 2][m + 2]

	/* verified against https://link.springer.com/content/pdf/bbm%3A978-1-4684-0208-7%2F1.pdf in rotation_tests.c*/

	Decimal cosp = cos(angle);
	Decimal sinp = sin(angle);
	Dw[0][0] = pow(cos(angle/2.), 4); // -2 -2
	Dw[1][0] = (-1/2.) * (1 + cosp) * sinp; // -1 -2
	Dw[2][0] = sqrt(3/8.) * pow(sinp, 2); // 0 -2
	Dw[3][0] = (1/2.) * sinp * (cosp - 1); // 1 -2
	Dw[4][0] = pow(sin(angle/2.), 4); // 2 -2

	Dw[0][1] = (1/2.) * (1 + cosp) * sinp; // -2 -1
	Dw[1][1] = pow(cosp, 2.) - (1/2.) * (1 - cosp); // -1 -1
	Dw[2][1] = -sqrt(3/8.) * sin(angle * 2.); // 0 -1
	Dw[3][1] = (1/2.) * (1 + cosp) - pow(cosp, 2); // 1 -1
	Dw[4][1] = (1/2.) * (cosp - 1) * sinp; // 2 -1

	Dw[0][2] = sqrt(3/8.) * pow(sinp, 2); // -2 0
	Dw[1][2] = sqrt(3/2.) * sinp * cosp; // -1 0
	Dw[2][2] = (1/2.) * (3 * pow(cosp, 2) - 1); // 0 0
	Dw[3][2] = -sqrt(3/2.) * sinp * cosp; // 1 0
	Dw[4][2] = sqrt(3/8.) * pow(sinp, 2);  // 2 0

	Dw[0][3] = -(1/2.) * (cosp - 1) * sinp; // -2 1 
	Dw[1][3] = (1/2.) * (1 + cosp) - pow(cosp, 2); // -1 1
	Dw[2][3] = sqrt(3/8.) * sin(angle * 2.); // 0 1
	Dw[3][3] = pow(cosp, 2) - (1/2.)*(1-cosp); // 1 1
	Dw[4][3] = (-1/2.) * (1 + cosp) * sinp; // 2 1

	Dw[0][4] = pow(sin(angle/2.), 4); // -2 2
	Dw[1][4] = (-1/2.) * (cosp - 1) * sinp; // -1 2
	Dw[2][4] = sqrt(3/8.) * pow(sinp, 2); // 0 2
	Dw[3][4] = (1/2.) * sinp * (cosp + 1); // 1 2  m' = 1, m = 2
	Dw[4][4] = pow(cos(angle / 2.), 4); // 2 2
}

/**
 * Function to free all allocated memory within Model.
 * @params m
 *  Pointer to Model to be freed.
 */
void free_all(struct Model *m) {
	unsigned int res, k;
	unsigned int params = m->params;

	for (res = 0; res < m->n_residues; res++) {
		if (m->residues[res].relaxation != NULL) {
			free(m->residues[res].relaxation);
			m->residues[res].relaxation = NULL;
		}
		if (m->residues[res].parameters != NULL) {
			free(m->residues[res].parameters);
			m->residues[res].parameters = NULL;
		}
		if (m->error_mode == 1) {
			if (m->residues[res].error_params != NULL) {
				for (k = 0; k < params; k++) {
					if (m->residues[res].error_params[k] != NULL) {
						free(m->residues[res].error_params[k]);
						m->residues[res].error_params[k] = NULL;
					}
				}
				free(m->residues[res].error_params);
				m->residues[res].error_params = NULL;
			}
			if (m->residues[res].errors_mean != NULL) {
				free(m->residues[res].errors_mean);
				m->residues[res].errors_mean = NULL;
			}
			if (m->residues[res].errors_std != NULL) {
				free(m->residues[res].errors_std);
				m->residues[res].errors_std = NULL;
			}
		}
	}
	free(m->residues);
	m->residues = NULL;

}

/**
 * Function to rotate spherical harmonics using Wigner D Matrices. I believe the
 * angles alpha, beta, gamma are aligned with the alpha, beta and gamma axis for deflections
 * (beta definitely is, I'm unsure about the others - but I'd assume so?)
 * @params or
 *  Pointer to Orient struct to modify.
 * @params alpha
 *  Alpha rotation
 * @params beta
 *  Beta rotation
 * @params gamma
 *  Gamma rotation
 * @returns Nothing
 *  Modifies orient inline.
 */
void rotate_Y2(struct Orient * or, Decimal alpha, Decimal beta, Decimal gamma) {
	// the idea should be to run calculate_Y2 _first_, and then this function.
	// it will take the orient and overwrite it with the rotated vectors.

	Decimal Dw[5][5];

	initialise_dwig(beta, Dw); // we initialise Dw[][] with the reduced Wigner matrices.


	
	/* the Wigner D-matrix D^j_{m', m} is e^(-i m' alpha) d^j_{m', m} e^(-i m gamma) */

	Complex Y2b[5]; // we store Y2.
	int i;
	Complex alpha_m[5]; // precalculation of the exp(-I * mp * alpha) and exp(-I * m * gamma) parts.
	Complex gamma_m[5];
	for (i = 0; i < 5; i++) {
		Y2b[i] = or->Y2[i];
		or->Y2[i] = 0;
		or->Y2c[i] = 0;
		alpha_m[i] = cexp(-(i - 2) * alpha * I);
		gamma_m[i] = cexp(-(i - 2) * gamma * I);
	}

	int m, mp;
	Complex Dcomp;

	for (mp = -2; mp <= 2; mp++) {
		for (m = -2; m <= 2; m++) {
			Dcomp = gamma_m[m+2];
			Dcomp *= Dw[mp+2][m+2];
			Dcomp *= alpha_m[mp+2]; 
			// D^j_{mp, m} = exp(-I * mp * alpha) * d^j_{m'm} * exp(-I * m * gamma)
			or->Y2[mp+2] += conj(Dcomp) * Y2b[m+2];
		}
		or->Y2c[mp+2] = conj(or->Y2[mp+2]); // and fill conjugate
	}

}

/**
 * Generates uniform random Decimal from 0 to 1 inclusive.
 * @return Decimal
 *  Long Decimal containing uniform random number.
 */
Decimal uniform_rand(void) {
	return ((Decimal) rand() + 1.) / ((Decimal) RAND_MAX + 1.);
}

void gen_params(const Decimal *minv, const Decimal *maxv, Decimal *pars, unsigned int n_pars) {
	unsigned int k;
	for (k = 0; k < n_pars; k++)
		pars[k] = minv[k] + (uniform_rand() * (maxv[k] - minv[k]));
}

void setup_paramlims(struct Model *m, Decimal S2NH, Decimal * minv, Decimal * maxv) {
	unsigned int k;
	for (k = 0; k < m->params; k++)
		minv[k] = 0;

	Decimal low_taus = 0.01;
	Decimal upp_taus = 100;
	Decimal low_tauf = 0;
	Decimal upp_tauf = 1;

	if (m->microsecond == ENABLED) {
		low_taus = 1000;
		upp_taus = 100000;
		low_tauf = 0;
		upp_tauf = 100;
	}
	switch (m->model) {
		case MOD_SMF:
			maxv[0] = upp_tauf; maxv[1] = 1;
			break;
		case MOD_SMFT:
			maxv[0] = upp_tauf; maxv[1] = 1; maxv[2] = 60000;
			break;
		case MOD_EMF:
			minv[0] = low_taus; minv[1] = S2NH;
			maxv[0] = upp_taus; maxv[1] = 1; maxv[2] = upp_tauf;
			break;
		case MOD_EMFT:
			minv[0] = low_taus; minv[1] = S2NH;
			maxv[0] = upp_taus; maxv[1] = 1; maxv[2] = upp_tauf; maxv[3] = 60000, maxv[4] = 60000;
			break;
		case MOD_DEMF:
			minv[0] = low_taus; minv[1] = S2NH; minv[3] = S2NH;
			maxv[0] = upp_taus; maxv[1] = 1; maxv[2] = upp_tauf; maxv[3] = 1;
			break;
		case MOD_DEMFT:
			minv[0] = low_taus; minv[1] = S2NH; minv[3] = S2NH;
			maxv[0] = upp_taus; maxv[1] = 1; maxv[2] = upp_tauf; maxv[3] = 1; maxv[4] = 60000, maxv[5] = 60000;
			break;
		case MOD_GAF:
			minv[0] = low_taus;
			maxv[0] = upp_taus; maxv[1] = upp_tauf;
			for (k = 2; k <= 7; k++) maxv[k] = 0.25;
			break;
		case MOD_GAFT:
			minv[0] = low_taus;
			maxv[0] = upp_taus; maxv[1] = upp_tauf;
			for (k = 2; k <= 7; k++) maxv[k] = 0.25;
			maxv[8] = 60000; maxv[9] = 60000;
			break;
		case MOD_AIMF:
			minv[0] = low_taus;
			maxv[0] = upp_taus; maxv[1] = upp_tauf;
			for (k = 2; k <= 7; k++) { maxv[k] = 1; minv[k] = S2NH; }
			break;
		case MOD_AIMFT:
			minv[0] = low_taus;
			maxv[0] = upp_taus; maxv[1] = upp_tauf;
			for (k = 2; k <= 7; k++) { maxv[k] = 1; minv[k] = S2NH; }
			maxv[8] = 60000; maxv[9] = 60000;
			break;
	    case MOD_EAIMF:
	        minv[0] = low_taus;
	        maxv[0] = upp_taus; maxv[1] = upp_tauf;
	        for (k = 2; k <= 5; k++) {maxv[k] = 1; minv[k] = S2NH; }
	        break;
	    case MOD_EAIMFT:
            minv[0] = low_taus;
            maxv[0] = upp_taus; maxv[1] = upp_tauf;
            for (k = 2; k <= 5; k++) {maxv[k] = 1; minv[k] = S2NH; }
            maxv[6] = 60000; maxv[7] = 60000;
            break;
	    case MOD_BGAF:
	        minv[0] = low_taus;
	        for (k = 1; k <= 3; k++) { maxv[k] = 0.25; }
	        break;
	    case MOD_BGAFT:
	        minv[0] = low_taus;
            for (k = 1; k <= 3; k++) { maxv[k] = 0.25; }
            maxv[4] = 60000;
            maxv[5] = 60000;
            break;
		case MOD_EGAF:
			minv[0] = low_taus;
			maxv[0] = upp_taus; maxv[1] = upp_tauf;
			for (k = 2; k <= 4; k++) maxv[k] = 0.25;
			minv[5] = S2NH; maxv[5] = 1;
			break;
		case MOD_EGAFT:
			minv[0] = low_taus;
			maxv[0] = upp_taus; maxv[1] = upp_tauf;
			for (k = 2; k <= 4; k++) maxv[k] = 0.25;
			minv[5] = S2NH; maxv[5] = 1;
			maxv[6] = 60000; maxv[7] = 60000;
			break;
		case MOD_BGF:
			minv[0] = low_taus;
			maxv[0] = upp_taus;
			maxv[1] = upp_tauf;
			break;
		case MOD_BGFT:
			minv[0] = low_taus;
			maxv[0] = upp_taus;
			maxv[1] = upp_tauf;
			maxv[8] = 60000;
			maxv[9] = 60000;
			break;
		default: break;
	}

	if (m->or_variation == VARIANT_A && m->ultrafast == ENABLED) {
		minv[m->params - 3] = 0; // alpha
		minv[m->params - 2] = 0; // beta
		minv[m->params - 1] = 0; // gamma
		maxv[m->params - 3] = 0; // alpha
		maxv[m->params - 2] = 0; // beta
		maxv[m->params - 1] = 0; // gamma
		minv[m->params - 4] = S2NH; // S2uf
		maxv[m->params - 4] = 1;
	} else if (m->or_variation == VARIANT_A) {
		minv[m->params - 3] = 0; // alpha
		minv[m->params - 2] = 0; // beta
		minv[m->params - 1] = 0; // gamma
		maxv[m->params - 3] = 0; // alpha
		maxv[m->params - 2] = 0; // beta
		maxv[m->params - 1] = 0; // gamma
	} else if (m->ultrafast == ENABLED) {
		minv[m->params - 1] = S2NH;
		maxv[m->params - 1] = 1;
	}
/**
		 * SMF parameters are \n
		 *   [0] tau\n
		 *   [1] S2\n
		 * SMFT parameters;\n
		 *   [0] tau\n
		 *   [1] S2\n
		 *   [2] Ea\n
		 * EMF parameters\n
		 *   [0] tau slow\n
		 *   [1] S2 slow\n
		 *   [2] tau fast\n
		 *   NOTE: The fast order parameter is calculated as S2NH/S2s\n
		 * EMFT parameters\n
		 *   [0] tau slow\n
		 *   [1] S2 slow\n
		 *   [2] tau fast\n
		 *   [3] activation energy for slow motion\n
		 *   [4] activation energy for fast motion\n
		 * DEMF parameters\n
		 *   [0] tau slow\n
		 *   [1] S2 slow\n
		 *   [2] tau fast\n
		 *   [3] S2 fast\n
		 * DEMFT parameters\n
		 *   [0] tau slow\n
		 *   [1] S2 slow\n
		 *   [2] tau fast\n
		 *   [3] S2 fast\n
		 *   [4] activation energy for slow motion\n
		 *   [5] activation energy for fast motion\n
		 * GAF parameters\n
		 *   [0] tau slow\n
		 *   [1] tau fast\n
		 *   [2-4] alpha, beta, gamma deflections for slow motions\n
		 *   [5-7] alpha, beta, gamma deflections for fast motions\n
		 * GAFT parameters\n
		 *   [0] tau slow\n
		 *   [1] tau fast\n
		 *   [2-4] alpha, beta, gamma deflections for slow motions\n
		 *   [5-7] alpha, beta, gamma deflections for fast motions\n
		 *   [8] activation energy for slow motion\n
		 *   [9] activation energy for fast motion\n
		 * EGAF parameters\n
		 *   [0] tau slow\n
		 *   [1] tau fast\n
		 *   [2-4] alpha, beta, gamma deflections for slow motions\n
		 *   [5] fast motion order parameter\n
		 * EGAFT parameters\n
		 *   [0] tau slow\n
		 *   [1] tau fast\n
		 *   [2-4] alpha, beta, gamma deflections for slow motions\n
		 *   [5] order parameter for fast motions\n
		 *   [6] activation energy for slow motion\n
		 *   [7] activation energy for fast motion\n
		 */
}

int determine_residues(unsigned int n_res, int myid, int numprocs, unsigned int *start, unsigned int *end) {
    // we have m->n_residues split over numprocs.
    // so


    unsigned int res_per_proc = (n_res / numprocs);
    if (n_res % numprocs != 0)
        res_per_proc++;
    *start = res_per_proc * myid;
    *end = *start + res_per_proc;
    if (*end > n_res)
        *end = n_res;
    return 1;
}
