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
double back_calc(long double * opts, struct Residue * resid, struct Relaxation * relax, struct Model * m, int *violations) {
	double calc_R = 0;
	unsigned int i;

	unsigned int model = m->model;
	if (relax->R <= 0)
		return -1;
	
	long double upper_lim_tf = (long double) 1 * pow(10, -8 + T_S);
	long double upper_lim_ts = (long double) 1 * pow(10, -5 + T_S);
	long double lower_lim_tf = (long double) 1 * pow(10, -13 + T_S);
	
	if (model == MOD_SMF || model == MOD_SMFT) {
		long double tau = opts[0];
		long double S2 = opts[1];
		long double tau_eff = 0, Ea = 0;
		
		if (tau < 0)
			(*violations)++;
		if (S2 < 0 || S2 > 1)
			(*violations)++;
		
		if (model == MOD_SMFT) {
			Ea = opts[2];
			tau_eff = tau * expl(Ea / (RYD * relax->T)); 
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
		long double taus = opts[0];
		long double S2s = opts[1];
		long double tauf = opts[2];
		long double taus_eff=0, tauf_eff=0, Eas=0, Eaf=0, S2f=resid->S2NH / S2s;

		if (model == MOD_EMFT) {
			Eas = opts[3];
			Eaf = opts[4];
		} else if (model == MOD_DEMFT) {
			Eas = opts[4];
			Eaf = opts[5];
		}
		if (model == MOD_EMFT || model == MOD_DEMFT) {
			taus_eff = taus * expl(Eas / (RYD * relax->T));
			tauf_eff = tauf * expl(Eaf / (RYD * relax->T));
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
		long double taus = opts[0];
		long double tauf = opts[1];
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		long double sigf[3] = {opts[5], opts[6], opts[7]};
		long double taus_eff=taus, tauf_eff=tauf, Eas=0, Eaf=0;
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
		long double taus = opts[0];
		long double tauf = opts[1];
		long double Ss[3] = {opts[2], opts[3], opts[4]};
		long double Sf[3] = {opts[5], opts[6], opts[7]};
		long double taus_eff=taus, tauf_eff=tauf, Eas=0, Eaf=0;
		/* WARNING: Need to perform variant orientations before entering back_calc.*/
		if (model == MOD_AIMFT) {
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
		long double taus = opts[0];
		long double tauf = opts[1];
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		long double S2f = opts[5];
		long double taus_eff = taus, tauf_eff = tauf, Eas=0, Eaf=0;
		if (model == MOD_EGAFT) {
			Eas = opts[6];
			Eaf = opts[7];
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
	long double papbN=0,papbC=0,kex=0, RDC=0, fm = 1;
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
		calc_R += (double) RDC;
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
double optimize_chisq(long double * opts, struct Residue * resid, struct Model * m, unsigned int params) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */

	unsigned int model = m->model;
	unsigned int or_variations = m->or_variation;
	double calc_R = 0;
	long double chisq = 0;
	int violations = 0;
	unsigned int i;
	double alpha, beta, gamma;
	if (or_variations == VARIANT_A) {
		alpha = (double) opts[params-3];
		beta = (double) opts[params-2];
		gamma = (double) opts[params-1];
		if (fabs(alpha) > OR_LIMIT || fabs(beta) > OR_LIMIT || fabs(gamma) > OR_LIMIT)
			violations++;
		for (i = 0; i < N_OR; i++) {
			calculate_Y2(&(resid->orients[i]));
			rotate_Y2(&(resid->orients[i]), alpha, beta, gamma);
		}
	}
	
	float mult = 1.;
	omp_set_num_threads(m->nthreads);
	#pragma omp parallel for reduction(+:chisq)
	for (i = 0; i < resid->n_relaxation; i++) {
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
	
	if (model == MOD_DEMF || model == MOD_DEMFT) {
		long double S2s = opts[1];
		long double S2f = opts[3];
		long double S2NH = S2s * S2f;
		
		chisq += ((pow(resid->S2NH - (double) S2NH, 2.)) / pow(resid->S2NHe, 2.));
	} else if (model == MOD_GAF || model == MOD_GAFT) {
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		long double sigf[3] = {opts[5], opts[6], opts[7]};
		
		double S2NHs, S2NHf, S2CHs, S2CHf, S2CNs, S2CNf; // S2CCs, S2CCf
		/** I'm unsure if the CC here is forward or backward so for now I have ignored it. */
		struct Orient *Os[] = {&(resid->orients[OR_NH]), &(resid->orients[OR_CNH]), &(resid->orients[OR_CN])};
		double *S2s[] = {&S2NHs, &S2CHs, &S2CNs};
		double *S2f[] = {&S2NHf, &S2CHf, &S2CNf};
		GAF_S2(sigs, Os, Os, S2s, 3, MODE_REAL);
		GAF_S2(sigf, Os, Os, S2f, 3, MODE_REAL);
		/* 100 weighting for order parameters */
		chisq += ((pow(resid->S2NH - (S2NHs * S2NHf), 2)) / pow(resid->S2NHe, 2));
		chisq += ((pow(resid->S2CH - (S2CHs * S2CHf), 2)) / pow(resid->S2CHe, 2));
		chisq += ((pow(resid->S2CN - (S2CNs * S2CNf), 2)) / pow(resid->S2CNe, 2));
	} else if (model == MOD_AIMF || model == MOD_AIMFT) {
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		long double sigf[3] = {opts[5], opts[6], opts[7]};
		
		double S2NHs, S2NHf, S2CHs, S2CHf, S2CNs, S2CNf; // S2CCs, S2CCf
		/** I'm unsure if the CC here is forward or backward so for now I have ignored it. */
		struct Orient *Os[] = {&(resid->orients[OR_NH]), &(resid->orients[OR_CNH]), &(resid->orients[OR_CN])};
		double *S2s[] = {&S2NHs, &S2CHs, &S2CNs};
		double *S2f[] = {&S2NHf, &S2CHf, &S2CNf};
		AIMF_S2(sigs, Os, S2s, 3);
		AIMF_S2(sigf, Os, S2f, 3);
		/* 100 weighting for order parameters */
		chisq += ((pow(resid->S2NH - (S2NHs * S2NHf), 2)) / pow(resid->S2NHe, 2));
		chisq += ((pow(resid->S2CH - (S2CHs * S2CHf), 2)) / pow(resid->S2CHe, 2));
		chisq += ((pow(resid->S2CN - (S2CNs * S2CNf), 2)) / pow(resid->S2CNe, 2));
	} else if (model == MOD_EGAF || model == MOD_EGAFT) {
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		double S2f = (double) opts[5];
		
		double S2NHs, S2CHs, S2CNs; // S2CCs
		/** I'm unsure if the CC here is forward or backward so for now I have ignored it. */
		struct Orient *Os[] = {&(resid->orients[OR_NH]), &(resid->orients[OR_CNH]), &(resid->orients[OR_CN])};
		double *S2s[] = {&S2NHs, &S2CHs, &S2CNs};
		GAF_S2(sigs, Os, Os, S2s, 3, MODE_REAL);
		/* 100 weghting for order parameters */
		chisq += ((pow(resid->S2NH - (S2NHs * S2f), 2)) / pow(resid->S2NHe, 2));
	//	chisq += 10*((pow(resid->S2CH - (S2CHs * S2f), 2)) / pow(resid->S2CHe, 2));
	//	chisq += 10*((pow(resid->S2CN - (S2CNs * S2f), 2)) / pow(resid->S2CNe, 2));
	}
	
	return (double) (chisq / resid->n_relaxation); 
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
int back_calculate(long double * opts, struct Residue * resid, struct Model * m, char *filename, unsigned int params) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */
	//if (opts[0] == -1)
	//	return -1;
	unsigned int or_variations = m->or_variation;
	double calc_R = 0;
	FILE * fp;
	fp = fopen(filename, "w");
	if (fp == NULL) {
		ERROR("%s not found.", filename);
		return -1;
	}
	int violations = 0;
	unsigned int i;
	double alpha, beta, gamma;
	if (or_variations == VARIANT_A) {
		alpha = (double) opts[params-3];
		beta = (double) opts[params-2];
		gamma = (double) opts[params-1];
		if (fabs(alpha) > OR_LIMIT || fabs(beta) > OR_LIMIT || fabs(gamma) > OR_LIMIT)
			violations++;
		for (i = 0; i < N_OR; i++) {
			calculate_Y2(&(resid->orients[i]));
			rotate_Y2(&(resid->orients[i]), alpha, beta, gamma);
		}
	}
	
	for (i = 0; i < resid->n_relaxation; i++) {
		calc_R = back_calc(opts, resid, &(resid->relaxation[i]), m, &violations);

		fprintf(fp, "%d\t%f\t%f\t%f", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);

		if (VERBOSE)
			fprintf(fp, "\t%f\t%f\t%f\t%f\t%d", resid->relaxation[i].field, resid->relaxation[i].wr, resid->relaxation[i].w1, resid->relaxation[i].T, resid->relaxation[i].type); 
		fprintf(fp, "\n");
	}


	fclose(fp);
	return -1;
}
