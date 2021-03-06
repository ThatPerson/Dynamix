/**
 * @file chisq.c
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>

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
Decimal back_calc(Decimal * opts, struct Residue * resid, struct Relaxation * relax, struct Model * m, int *violations) {
	Decimal calc_R = 0;
	unsigned int i;

	unsigned int model = m->model;
	if (relax->R <= 0)
		return -1;
	
	Decimal upper_lim_tf = (Decimal) 1 * pow(10, -8 + T_S);
	Decimal upper_lim_ts = (Decimal) 1 * pow(10, -5 + T_S);
	Decimal lower_lim_tf = (Decimal) 1 * pow(10, -13 + T_S);
	
	if (model == MOD_SMF || model == MOD_SMFT) {
		Decimal tau = opts[0];
		Decimal S2 = opts[1];
		Decimal tau_eff = 0, Ea = 0;
		
		if (tau < 0)
			(*violations)++;
		if (S2 < 0 || S2 > 1)
			(*violations)++;
		
		if (model == MOD_SMFT) {
			Ea = opts[2];
			tau_eff = tau * exp(Ea / (RYD * relax->T));
		} else {
			tau_eff = tau;			
		}
		
		switch (relax->type) {
			case R_15NR1:
				calc_R = SMF_R1(resid, relax, tau_eff, S2, MODE_15N);
				break;
			case R_15NR1p:
				calc_R = SMF_R2(resid, relax, tau_eff, S2, MODE_15N);
				break;
			case R_13CR1:
				calc_R = SMF_R1(resid, relax, tau_eff, S2, MODE_13C);
				break;
			case R_13CR1p:
				calc_R = SMF_R2(resid, relax, tau_eff, S2, MODE_13C);
				break;
			default:
				ERROR("Unknown relaxation type: %d", relax->type);
				exit(0);
				break;
		}
	} else if (model == MOD_EMF || model == MOD_EMFT || model == MOD_DEMF || model == MOD_DEMFT) {
		Decimal taus = opts[0];
		Decimal S2s = opts[1];
		Decimal tauf = opts[2];
		Decimal taus_eff=0, tauf_eff=0, Eas=0, Eaf=0, S2f=resid->S2NH / S2s;

		if (model == MOD_EMFT) {
			Eas = opts[3];
			Eaf = opts[4];
		} else if (model == MOD_DEMFT) {
			Eas = opts[4];
			Eaf = opts[5];
		}
		if (model == MOD_EMFT || model == MOD_DEMFT) {
			taus_eff = taus * exp(Eas / (RYD * relax->T));
			tauf_eff = tauf * exp(Eaf / (RYD * relax->T));
		} else {
			taus_eff = taus;
			tauf_eff = tauf;
		}
		if (model == MOD_DEMF || model == MOD_DEMFT)
			S2f = opts[3];
		
		if (taus < 0 || tauf < 0)
			(*violations)++;
		if (S2s < resid->S2NH || S2f < resid->S2NH || S2s > 1 || S2f > 1)
			(*violations)++;
		if (tauf_eff > taus_eff)
			(*violations)++;
		if (tauf_eff > upper_lim_tf || taus_eff > upper_lim_ts || tauf_eff < lower_lim_tf)
			(*violations)++;
		
		switch (relax->type) {
			case R_15NR1:
				calc_R = EMF_R1(resid, relax, taus_eff, S2s, tauf_eff, S2f, MODE_15N);
				break;
			case R_15NR1p:
				calc_R = EMF_R2(resid, relax, taus_eff, S2s, tauf_eff, S2f, MODE_15N);
				break;
			case R_13CR1:
				calc_R = EMF_R1(resid, relax, taus_eff, S2s, tauf_eff, S2f, MODE_13C);
				break;
			case R_13CR1p:
				calc_R = EMF_R2(resid, relax, taus_eff, S2s, tauf_eff, S2f, MODE_13C);
				break;
			default:
				ERROR("Unknown relaxation type: %d", relax->type);
				exit(0);
				break;
		}
	} else if (model == MOD_GAF || model == MOD_GAFT) {
		Decimal taus = opts[0];
		Decimal tauf = opts[1];
		Decimal sigs[3] = {opts[2], opts[3], opts[4]};
		Decimal sigf[3] = {opts[5], opts[6], opts[7]};
		Decimal taus_eff=taus, tauf_eff=tauf, Eas=0, Eaf=0;
		/* WARNING: Need to perform variant orientations before entering back_calc.*/
		if (model == MOD_GAFT) {
			Eas = opts[8];
			Eaf = opts[9];
			taus_eff *= expl(Eas / (RYD * relax->T));
			tauf_eff *= expl(Eaf / (RYD * relax->T));
		}
		
		if (taus < 0 || tauf < 0)
			(*violations)++;
		if (tauf_eff > taus_eff)
			(*violations)++;
		if (tauf_eff > upper_lim_tf || taus_eff > upper_lim_ts || tauf_eff < lower_lim_tf)
			(*violations)++;
		for (i = 0; i < 3; i++) {
			if (sigs[i] < 0 || sigs[i] > 0.52360)
				(*violations)++;
			if (sigf[i] < 0 || sigf[i] > 0.52360)
				(*violations)++;
		}
		
		switch (relax->type) {
			case R_15NR1:
				calc_R = GAF_15NR1(resid, relax, taus_eff, tauf_eff, sigs, sigf);
				break;
			case R_15NR1p:
				calc_R = GAF_15NR2(resid, relax, taus_eff, tauf_eff, sigs, sigf);
				break;
			case R_13CR1:
				calc_R = GAF_13CR1(resid, relax, taus_eff, tauf_eff, sigs, sigf);
				break;
			case R_13CR1p:
				calc_R = GAF_13CR2(resid, relax, taus_eff, tauf_eff, sigs, sigf);
				break;
			default:
				ERROR("Unknown relaxation type: %d", relax->type);
				exit(0);
				break;
		}
	} else if (model == MOD_AIMF || model == MOD_AIMFT) {
		Decimal taus = opts[0];
		Decimal tauf = opts[1];
		Decimal Ss[3] = {opts[2], opts[3], opts[4]};
		Decimal Sf[3] = {opts[5], opts[6], opts[7]};
		Decimal taus_eff=taus, tauf_eff=tauf, Eas=0, Eaf=0;
		/* WARNING: Need to perform variant orientations before entering back_calc.*/
		if (model == MOD_AIMFT) {
			Eas = opts[8];
			Eaf = opts[9];
			taus_eff *= exp(Eas / (RYD * relax->T));
			tauf_eff *= exp(Eaf / (RYD * relax->T));
		}
		
		if (taus < 0 || tauf < 0)
			(*violations)++;
		if (tauf_eff > taus_eff)
			(*violations)++;
		if (tauf_eff > upper_lim_tf || taus_eff > upper_lim_ts || tauf_eff < lower_lim_tf)
			(*violations)++;
		for (i = 0; i < 3; i++) {
			if (Ss[i] < 0 || Ss[i] > 1)
				(*violations)++;
			if (Sf[i] < 0 || Sf[i] > 1)
				(*violations)++;
		}
		
		switch (relax->type) {
			case R_15NR1:
				calc_R = AIMF_15NR1(resid, relax, taus_eff, tauf_eff, Ss, Sf);
				break;
			case R_15NR1p:
				calc_R = AIMF_15NR2(resid, relax, taus_eff, tauf_eff, Ss, Sf);
				break;
			case R_13CR1:
				calc_R = AIMF_13CR1(resid, relax, taus_eff, tauf_eff, Ss, Sf);
				break;
			case R_13CR1p:
				calc_R = AIMF_13CR2(resid, relax, taus_eff, tauf_eff, Ss, Sf);
				break;
			default:
				ERROR("Unknown relaxation type: %d", relax->type);
				exit(0);
				break;
		}
	} else if (model == MOD_EGAF || model == MOD_EGAFT) {
		Decimal taus = opts[0];
		Decimal tauf = opts[1];
		Decimal sigs[3] = {opts[2], opts[3], opts[4]};
		Decimal S2f = opts[5];
		Decimal taus_eff = taus, tauf_eff = tauf, Eas=0, Eaf=0;
		if (model == MOD_EGAFT) {
			Eas = opts[6];
			Eaf = opts[7];
			taus_eff *= exp(Eas / (RYD * relax->T));
			tauf_eff *= exp(Eaf / (RYD * relax->T));
		}
		
		if (taus < 0 || tauf < 0)
			(*violations)++;
		if (tauf_eff > taus_eff)
			(*violations)++;
		if (tauf_eff > upper_lim_tf || taus_eff > upper_lim_ts || tauf_eff < lower_lim_tf)
			(*violations)++;
		for (i = 0; i < 3; i++) {
			if (sigs[i] < 0 || sigs[i] > 0.52360)
				(*violations)++;
		}
		if (S2f < resid->S2NH || S2f > 1)
			(*violations)++;
		
		switch (relax->type) {
			case R_15NR1:
				calc_R = EGAF_15NR1(resid, relax, taus_eff, tauf_eff, sigs, S2f);
				break;
			case R_15NR1p:
				calc_R = EGAF_15NR2(resid, relax, taus_eff, tauf_eff, sigs, S2f);
				break;
			case R_13CR1:
				calc_R = EGAF_13CR1(resid, relax, taus_eff, tauf_eff, sigs, S2f);
				break;
			case R_13CR1p:
				calc_R = EGAF_13CR2(resid, relax, taus_eff, tauf_eff, sigs, S2f);
				break;
			default:
				ERROR("Unknown relaxation type: %d", relax->type);
				exit(0);
				break;
		}
	} else {
		printf("Model not yet implemented.\n");
		calc_R = -1;
	}
	Decimal papbN=0,papbC=0,kex=0, RDC=0, fm;
	if (m->rdc == RDC_ON) {
		if (m->or_variation == VARIANT_A) {
			papbC = opts[m->params - 6];
			papbN = opts[m->params - 5];
			kex  = opts[m->params - 4];
		} else {
			papbC = opts[m->params - 3];
			papbN = opts[m->params - 2];
			kex  = opts[m->params - 1];
		}

		if (papbC < 0)
			(*violations)++;
		if (papbN < 0)
			(*violations)++;
		if (kex < 0)
			(*violations)++;
		
		fm = relax->field / 700.;
		if (relax->type == R_13CR1p) {
			fm *= 9.869683408806043/3.976489314034722;
			RDC = papbC * sq(fm) * kex;
			RDC /= (sq(relax->w1) + sq(kex));
		} else if (relax->type == R_15NR1p) {
			fm *= 1;
			RDC = papbN * sq(fm) * kex;
			RDC /= (sq(relax->w1) + sq(kex));
		} else {
			RDC = 0;
		}
		calc_R += (Decimal) RDC;
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
	Decimal calc_R = 0;
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
	omp_set_num_threads(m->nthreads);
	#pragma omp parallel for reduction(+:chisq)
	for (i = 0; i < resid->n_relaxation; i++) {
		if (resid->relaxation[i].R == -1)
			continue;
		calc_R = back_calc(opts, resid, &(resid->relaxation[i]), m, &violations);
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
	
	Decimal s2_mult = 100;
	
	if (model == MOD_DEMF || model == MOD_DEMFT) {
		Decimal S2s = opts[1];
		Decimal S2f = opts[3];
		Decimal S2NH = S2s * S2f;
		
		chisq += m->WS2NH * ((pow(resid->S2NH - (Decimal) S2NH, 2.)) / pow(resid->S2NHe, 2.));
	} else if (model == MOD_GAF || model == MOD_GAFT) {
		Decimal sigs[3] = {opts[2], opts[3], opts[4]};
		Decimal sigf[3] = {opts[5], opts[6], opts[7]};
		
		Decimal S2NHs, S2NHf, S2CHs, S2CHf, S2CNs, S2CNf, S2CCs, S2CCf;
		/** I'm unsure if the CC here is forward or backward so for now I have ignored it. */
		struct Orient *Os[] = {&(resid->orients[OR_NH]), &(resid->orients[OR_CNH]), &(resid->orients[OR_CN]), &(resid->orients[OR_CCAp])};
		// CCAp is shorter than CCAc.
		Decimal *S2s[] = {&S2NHs, &S2CHs, &S2CNs, &S2CCs};
		Decimal *S2f[] = {&S2NHf, &S2CHf, &S2CNf, &S2CCf};
		GAF_S2(sigs, Os, Os, S2s, 4, MODE_REAL);
		GAF_S2(sigf, Os, Os, S2f, 4, MODE_REAL);
		/* 100 weighting for order parameters */
		chisq += m->WS2NH * ((pow(resid->S2NH - (S2NHs * S2NHf), 2)) / pow(resid->S2NHe, 2));
		if (m->WS2CH != 0) chisq += m->WS2CH * ((pow(resid->S2CH - (S2CHs * S2CHf), 2)) / pow(resid->S2CHe, 2));
		if (m->WS2CN != 0) chisq += m->WS2CN * ((pow(resid->S2CN - (S2CNs * S2CNf), 2)) / pow(resid->S2CNe, 2));
		if (m->WS2CC != 0) chisq += m->WS2CC * ((pow(resid->S2CC - (S2CCs * S2CCf), 2)) / pow(resid->S2CCe, 2));
	} else if (model == MOD_AIMF || model == MOD_AIMFT) {
		Decimal sigs[3] = {opts[2], opts[3], opts[4]};
		Decimal sigf[3] = {opts[5], opts[6], opts[7]};
		
		Decimal S2NHs, S2NHf, S2CHs, S2CHf, S2CNs, S2CNf, S2CCs, S2CCf; // S2CCs, S2CCf
	
		struct Orient *Os[] = {&(resid->orients[OR_NH]), &(resid->orients[OR_CNH]), &(resid->orients[OR_CN]), &(resid->orients[OR_CCAp])};
		Decimal *S2s[] = {&S2NHs, &S2CHs, &S2CNs, &S2CCs};
		Decimal *S2f[] = {&S2NHf, &S2CHf, &S2CNf, &S2CCf};
		AIMF_S2(sigs, Os, S2s, 4);
		AIMF_S2(sigf, Os, S2f, 4);
		/* 100 weighting for order parameters */
		chisq += s2_mult * ((pow(resid->S2NH - (S2NHs * S2NHf), 2)) / pow(resid->S2NHe, 2));
		if (m->WS2CH != 0) chisq += m->WS2CH * ((pow(resid->S2CH - (S2CHs * S2CHf), 2)) / pow(resid->S2CHe, 2));
		if (m->WS2CN != 0) chisq += m->WS2CN * ((pow(resid->S2CN - (S2CNs * S2CNf), 2)) / pow(resid->S2CNe, 2));
		if (m->WS2CC != 0) chisq += m->WS2CC * ((pow(resid->S2CC - (S2CCs * S2CCf), 2)) / pow(resid->S2CCe, 2));
	} else if (model == MOD_EGAF || model == MOD_EGAFT) {
		Decimal sigs[3] = {opts[2], opts[3], opts[4]};
		Decimal S2f = (Decimal) opts[5];
		
		Decimal S2NHs, S2CHs, S2CNs, S2CCs; // S2CCs

		struct Orient *Os[] = {&(resid->orients[OR_NH]), &(resid->orients[OR_CNH]), &(resid->orients[OR_CN]), &(resid->orients[OR_CCAp])};
		Decimal *S2s[] = {&S2NHs, &S2CHs, &S2CNs, &S2CCs};
		GAF_S2(sigs, Os, Os, S2s, 4, MODE_REAL);
		/* 100 weghting for order parameters */
		chisq += m->WS2NH * ((pow(resid->S2NH - (S2NHs * S2f), 2)) / pow(resid->S2NHe, 2));
	//	chisq += 10*((pow(resid->S2CH - (S2CHs * S2f), 2)) / pow(resid->S2CHe, 2));
	//	chisq += 10*((pow(resid->S2CN - (S2CNs * S2f), 2)) / pow(resid->S2CNe, 2));
	//	chisq += 10*((pow(resid->S2CC - (S2CCs * S2f), 2)) / pow(resid->S2CCe, 2));
	}
	
	return (Decimal) (chisq / resid->n_relaxation); 
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
	Decimal calc_R = 0;
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
	
	for (i = 0; i < resid->n_relaxation; i++) {
		calc_R = back_calc(opts, resid, &(resid->relaxation[i]), m, &violations);

		fprintf(fp, "%d\t%lf\t%lf\t%lf", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);

		if (VERBOSE)
			fprintf(fp, "\t%lf\t%lf\t%lf\t%lf\t%d", resid->relaxation[i].field, resid->relaxation[i].wr, resid->relaxation[i].w1, resid->relaxation[i].T, resid->relaxation[i].type); 
		fprintf(fp, "\n");
	}


	fclose(fp);
	return -1;
}
