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
 * chisq.c
 *
 * Provides functions to backcalculate relaxation rate parameters, and calculate a
 * reduced chisq value from these to allow for parameter optimization.
 */

#include <math.h>
#include <stdio.h>
#include "datatypes.h"
#include "models/model.h"

int bcpars_init(struct BCParameters *pars, Decimal slow, Decimal fast) {
    pars->taus = 0;
    pars->tauf = 0;
    pars->Eas = -1;
    pars->Eaf = -1;
    pars->S2NHs = slow;
    pars->S2NHf = fast;
    pars->S2NCSAxs = slow;
    pars->S2NCSAxf = fast;
    pars->S2NCSAys = slow;
    pars->S2NCSAyf = fast;
    pars->S2NCSAxys = slow;
    pars->S2NCSAxyf = fast;
    pars->S2CSAxs = slow;
    pars->S2CSAxf = fast;
    pars->S2CSAys = slow;
    pars->S2CSAyf = fast;
    pars->S2CSAxys = slow;
    pars->S2CSAxyf = fast;
    pars->S2CNs = slow;
    pars->S2CNf = fast;
    pars->S2CaNs = slow;
    pars->S2CaNf = fast;
    pars->S2CHs = slow;
    pars->S2CHf = fast;
    pars->S2CCAps = slow;
    pars->S2CCApf = fast;
    pars->S2CCAcs = slow;
    pars->S2CCAcf = fast;
    pars->S2NHrs = slow;
    pars->S2NHrf = fast;
    pars->S2CHrs = slow;
    pars->S2CHrf = fast;
    return 0;
}

void check_S2_violations(struct BCParameters *pars, int *violations) {
    if (pars->S2NHs > 1 || pars->S2NHs < 0)
        (*violations)++;
    if (pars->S2NHf > 1 || pars->S2NHf < 0)
        (*violations)++;
    if (pars->S2CCAps > 1 || pars->S2CCAps < 0)
        (*violations)++;
    if (pars->S2CCApf > 1 || pars->S2CCApf < 0)
        (*violations)++;
    if (pars->S2CHs > 1 || pars->S2CHs < 0)
        (*violations)++;
    if (pars->S2CHf > 1 || pars->S2CHf < 0)
        (*violations)++;
    if (pars->S2CNs > 1 || pars->S2CNs < 0)
        (*violations)++;
    if (pars->S2CNf > 1 || pars->S2CNf < 0)
        (*violations)++;
}

int opts_to_bcpars(Decimal *opts, struct BCParameters *pars, unsigned int model, struct Residue *resid, int *violations) {
    Decimal S2s, S2f;
    int i;
    Decimal *S2sP[] = {&(pars->S2NHs),
                      &(pars->S2NCSAxs),
                      &(pars->S2NCSAys),
                      &(pars->S2NCSAxys),
                      &(pars->S2CSAxs),
                      &(pars->S2CSAys),
                      &(pars->S2CSAxys),
                      &(pars->S2CNs),
                      &(pars->S2CaNs),
                      &(pars->S2CHs),
                      &(pars->S2CCAcs),
                      &(pars->S2CCAps)
    };
    Decimal *S2fP[] = {&(pars->S2NHf),
                      &(pars->S2NCSAxf),
                      &(pars->S2NCSAyf),
                      &(pars->S2NCSAxyf),
                      &(pars->S2CSAxf),
                      &(pars->S2CSAyf),
                      &(pars->S2CSAxyf),
                      &(pars->S2CNf),
                      &(pars->S2CaNf),
                      &(pars->S2CHf),
                      &(pars->S2CCAcf),
                      &(pars->S2CCApf)
    };
    struct Orient *As[] = {
            &(resid->orients[OR_NH]),
            &(resid->orients[OR_NCSAxx]),
            &(resid->orients[OR_NCSAyy]),
            &(resid->orients[OR_NCSAxx]),
            &(resid->orients[OR_CCSAxx]),
            &(resid->orients[OR_CCSAyy]),
            &(resid->orients[OR_CCSAxx]),
            &(resid->orients[OR_CN]),
            &(resid->orients[OR_NCA]),
            &(resid->orients[OR_CNH]),
            &(resid->orients[OR_CCAc]),
            &(resid->orients[OR_CCAp])
    };
    struct Orient *Bs[] = {
            &(resid->orients[OR_NH]),
            &(resid->orients[OR_NCSAxx]),
            &(resid->orients[OR_NCSAyy]),
            &(resid->orients[OR_NCSAyy]),
            &(resid->orients[OR_CCSAxx]),
            &(resid->orients[OR_CCSAyy]),
            &(resid->orients[OR_CCSAyy]),
            &(resid->orients[OR_CN]),
            &(resid->orients[OR_NCA]),
            &(resid->orients[OR_CNH]),
            &(resid->orients[OR_CCAc]),
            &(resid->orients[OR_CCAp])
    };
    if (model == MOD_SMF || model == MOD_SMFT) {
        // tau, S2, [Ea]
        bcpars_init(pars, 1, opts[1]);
        pars->tauf = opts[0];
        if (model == MOD_SMFT)
            pars->Eaf = opts[2];
    } else if (model == MOD_EMF || model == MOD_DEMF || model == MOD_EMFT || model == MOD_DEMFT){
        S2s = opts[1];
        S2f = resid->S2NH / S2s;
        if (model == MOD_DEMFT || model == MOD_DEMF) {
            S2f = opts[3];
        }
        bcpars_init(pars, S2s, S2f);

        pars->taus = opts[0];
        pars->tauf = opts[2];
        if (model == MOD_DEMFT) {
            pars->Eas = opts[4];
            pars->Eaf = opts[5];
        } else if (model == MOD_EMFT) {
            pars->Eas = opts[3];
            pars->Eaf = opts[4];
        }
    } else if (model == MOD_GAF || model == MOD_GAFT) {
        bcpars_init(pars, 0, 0);
        pars->taus = opts[0];
        pars->tauf = opts[1];
        Decimal sigs[3] = {opts[2], opts[3], opts[4]};
        Decimal sigf[3] = {opts[5], opts[6], opts[7]};
        for (i = 0; i < 3; i++) {
            if (sigs[i] < 0 || sigs[i] > 0.52360)
                (*violations)++;
            if (sigf[i] < 0 || sigf[i] > 0.52360)
                (*violations)++;
        }
        GAF_S2(sigs, As, Bs, S2sP, 12, MODE_REAL);
        GAF_S2(sigf, As, Bs, S2fP, 12, MODE_REAL);
        if (model == MOD_GAFT) {
            pars->Eas = opts[8];
            pars->Eaf = opts[9];
        }
        pars->S2NHrs = 1; pars->S2NHrf = 1; pars->S2CHrs = 1; pars->S2CHrf = 1;
    } else if (model == MOD_EGAF || model == MOD_EGAFT) {
        Decimal sigs[3] = {opts[2], opts[3], opts[4]};
        S2f = opts[5];
        for (i = 0; i < 3; i++) {
            if (sigs[i] < 0 || sigs[i] > 0.52360)
                (*violations)++;
        }
        bcpars_init(pars, 0, S2f);
        GAF_S2(sigs, As, Bs, S2sP, 12, MODE_REAL);
        pars->taus = opts[0];
        pars->tauf = opts[1];
        if (model == MOD_EGAFT) {
            pars->Eas = opts[6];
            pars->Eaf = opts[7];
        }
        pars->S2NHrs = 1; pars->S2NHrf = 1; pars->S2CHrs = 1; pars->S2CHrf = 1;
    } else if (model == MOD_AIMF || model == MOD_AIMFT) {
        bcpars_init(pars, 0, 0);
        pars->taus = opts[0];
        pars->tauf = opts[1];
        Decimal sigs[3] = {opts[2], opts[3], opts[4]};
        Decimal sigf[3] = {opts[5], opts[6], opts[7]};
        AIMF_S2(sigs, As, S2sP, 12);
        AIMF_S2(sigf, As, S2fP, 12);
        if (model == MOD_AIMFT) {
            pars->Eas = opts[8];
            pars->Eaf = opts[9];
        }
        pars->S2NHrs = 1; pars->S2NHrf = 1; pars->S2CHrs = 1; pars->S2CHrf = 1;
    } else {
        ERROR("Model %d does not exist.\n", model);
        return 1;
    }

    Decimal upper_lim_tf = (Decimal) 1 * pow(10, -8 + T_S);
    Decimal upper_lim_ts = (Decimal) 1 * pow(10, -5 + T_S);
    Decimal lower_lim_tf = (Decimal) 1 * pow(10, -13 + T_S);

    Decimal ttaus = pars->taus, ttauf = pars->tauf;
    if (pars->Eas != -1 || pars->Eaf != -1) {
        ttaus = temp_tau(pars->taus, pars->Eas, 300);
        ttauf = temp_tau(pars->tauf, pars->Eaf, 300);
    }

    check_S2_violations(pars, violations);

   // printf("Pars %lf (%lf)\n", ttaus, upper_lim_ts);
    if (ttaus < ttauf)
        (*violations)++;
    if (ttaus > upper_lim_ts)
        (*violations)++;
    if (ttauf > upper_lim_tf)
        (*violations)++;
    if (ttauf < lower_lim_tf)
        (*violations)++;
    for (i = 0; i < 11; i++) {
        if (*(S2sP[0]) > 1 || (*S2sP[0]) < 0)
            (*violations)++;
    }
    return 0;
}


/**
 * Relaxation back calculation. Takes parameters, residue and relaxation and returns
 * calculated relaxation rate.
 * @params *opts
 *  Array of parameters
 * @params *resid
 *  Pointer to residue
 * @params *relax
 *  Pointer to relaxation rate
 * @params model
 *  Model definition
 * @params *violations
 *  Int returning number of violations of constraints (for chisq)
 * @return calculated R value
 */
Decimal back_calc(struct Residue * resid, struct Relaxation * relax, struct Model * m, const int *violations, struct BCParameters *pars) {
	Decimal calc_R = -1;
    (void) violations;
    (void) m;
	if (relax->R <= 0)
		return -1;

    switch (relax->type) {
        case R_15NR1:  calc_R = Calc_15NR1(resid, relax, pars, m); break;
        case R_15NR1p: calc_R = Calc_15NR2(resid, relax, pars, m); break;
        case R_13CR1:  calc_R = Calc_13CR1(resid, relax, pars, m); break;
        case R_13CR1p: calc_R = Calc_13CR2(resid, relax, pars, m); break;
        default: ERROR("Relaxation type %d unknown\n", relax->type); break;
    }

	return calc_R;
}

/**
 * Optimization function. Loops over all relaxation measurements and predicts using models.\n
 * Calculates (real - calc)^2 / (error^2) over all relaxation measurements and returns.\n
 * @param opts
 *  Pointer to array containing parameters.
 * @param resid
 *  Pointer to residue being considered
 * @param model
 *  MOD_SMF etc.
 * @return Returns chisq value.
 */
Decimal optimize_chisq(Decimal * opts, struct Residue * resid, struct Model * m, unsigned int params) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */

	unsigned int model = m->model;
	unsigned int or_variations = m->or_variation;
	Decimal calc_R;
	Decimal chisq = 0;
	int violations = 0;
	unsigned int i;
	Decimal alpha, beta, gamma;
	if (or_variations == VARIANT_A) {
		alpha = (Decimal) opts[params-3];
		beta = (Decimal) opts[params-2];
		gamma = (Decimal) opts[params-1];
		if (fabs(alpha) > OR_LIMIT || fabs(beta) > OR_LIMIT || fabs(gamma) > OR_LIMIT)
			violations++;
		for (i = 0; i < N_OR; i++) {
			calculate_Y2(&(resid->orients[i]));
			rotate_Y2(&(resid->orients[i]), alpha, beta, gamma);
		}
	}
	
	Decimal mult = 1.;
	struct BCParameters pars;
    int k = opts_to_bcpars(opts, &pars, model, resid, &violations);
    if (k != 0)
        return -1;
#pragma omp parallel for reduction(+:chisq)
	for (i = 0; i < resid->n_relaxation; i++) {
		if (resid->relaxation[i].R == -1)
			continue;
		calc_R = back_calc(resid, &(resid->relaxation[i]), m, &violations, &pars);
		if (m->cn_ratio == CNRATIO_ON) {
			mult = 1;
			if (resid->relaxation[i].type == R_13CR1 && resid->relaxation[i].type == R_13CR1p)
				mult = resid->cn;
		}
		
		chisq += mult*((pow(resid->relaxation[i].R - calc_R, 2.)) / pow(resid->relaxation[i].Rerror, 2.));
	}
	//dp += resid->n_relaxation;
	
	chisq += 100000 * violations;
	if (violations > 0) {
		LOG("Chisq optimization had %d violations.", violations);
	}
	
    Decimal S2NH, S2CH, S2CN, S2CC;
	S2NH = pars.S2NHs * pars.S2NHf;
	S2CH = pars.S2CHs * pars.S2CHf;
	S2CN = pars.S2CNs * pars.S2CNf;
	S2CC = pars.S2CCAps * pars.S2CCApf;

    chisq += m->WS2NH * ((pow(resid->S2NH - S2NH, 2.)) / pow(resid->S2NHe, 2.));
    chisq += m->WS2CH * ((pow(resid->S2CH - S2CH, 2.)) / pow(resid->S2CHe, 2.));
    chisq += m->WS2CN * ((pow(resid->S2CN - S2CN, 2.)) / pow(resid->S2CNe, 2.));
    chisq += m->WS2CC * ((pow(resid->S2CC - S2CC, 2.)) / pow(resid->S2CCe, 2.));
    unsigned int div = resid->n_relaxation + m->WS2CC + m->WS2CH + m->WS2CN + m->WS2NH;

	return (Decimal) (chisq / div);
	/* normalise to number of relaxation measurements - otherwise when using like 85 the chisq becomes huge which hinders convergence */
}

/**
 * Back calculation function. Loops over all relaxation measurements and predicts using models.\n
 * Outputs file containing [relaxation id, calculated, real, real error].\n
 * @param opts
 *  Pointer to array containing parameters.
 * @param resid
 *  Pointer to residue being considered
 * @param model
 *  MOD_SMF etc.
 * @param filename
 *  File to output calculations into
 * @return Returns 1 if successful, else -1.
 */
int back_calculate(Decimal * opts, struct Residue * resid, struct Model * m, char *filename, unsigned int params) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */
	//if (opts[0] == -1)
	//	return -1;
	unsigned int or_variations = m->or_variation;
	Decimal calc_R;
	FILE * fp;
	fp = fopen(filename, "w");
	if (fp == NULL) {
		ERROR("%s not found.", filename);
		return -1;
	}
	int violations = 0;
	unsigned int i;
	Decimal alpha, beta, gamma;
	if (or_variations == VARIANT_A) {
		alpha = (Decimal) opts[params-3];
		beta = (Decimal) opts[params-2];
		gamma = (Decimal) opts[params-1];
		if (fabs(alpha) > OR_LIMIT || fabs(beta) > OR_LIMIT || fabs(gamma) > OR_LIMIT)
			violations++;
		for (i = 0; i < N_OR; i++) {
			calculate_Y2(&(resid->orients[i]));
			rotate_Y2(&(resid->orients[i]), alpha, beta, gamma);
		}
	}
    struct BCParameters pars;
    int k = opts_to_bcpars(opts, &pars, m->model, resid, &violations);
    if (k != 0) {
        return -1;
    }
	for (i = 0; i < resid->n_relaxation; i++) {
		calc_R = back_calc(resid, &(resid->relaxation[i]), m, &violations, &pars);

		fprintf(fp, "%d\t%lf\t%lf\t%lf", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);

		if (VERBOSE)
			fprintf(fp, "\t%lf\t%lf\t%lf\t%lf\t%d", resid->relaxation[i].field, resid->relaxation[i].wr, resid->relaxation[i].w1, resid->relaxation[i].T, resid->relaxation[i].type); 
		fprintf(fp, "\n");
	}


	fclose(fp);
	return -1;
}
