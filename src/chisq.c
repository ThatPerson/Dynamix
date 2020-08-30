/**
 * @file chisq.c
 */

#include <stdio.h>
#include <math.h>

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
double optimize_chisq(long double * opts, struct Residue * resid, int model) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */
	double calc_R = 0;
	double chisq = 0;
	long double upper_lim_tf = (long double) 1 * pow(10, -8 + T_S);
	long double upper_lim_ts = (long double) 1 * pow(10, -5 + T_S);
	long double lower_lim_tf = (long double) 1 * pow(10, -13 + T_S);
	//printf("%d %d: %d\n", model, MOD_SMFT, model == MOD_SMFT);
	if (model == MOD_SMF || model == MOD_SMFT) {
		/* Simple Model Free Analysis */
		long double tau = opts[0];
		long double S2 = opts[1];
		long double tau_eff=0;
		long double Ea=0;
		if (model == MOD_SMFT)
			Ea = opts[2];
		else
			tau_eff = tau;

		if (tau < 0)
			chisq += 100000000;
		if (S2 < 0 || S2 > 1)
			chisq += 100000000;
		int i;
		for (i = 0; i < resid->n_relaxation; i++) {
			if (resid->relaxation[i].R <= 0)
                                continue;

			if (model == MOD_SMFT)
				tau_eff = tau * expl(Ea / (RYD * resid->relaxation[i].T));
			//printf("%Le\n", tau_eff);
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = SMF_R1(resid, &(resid->relaxation[i]), tau_eff, S2, MODE_15N);
					break;
				case R_15NR1p:
					calc_R = SMF_R2(resid, &(resid->relaxation[i]), tau_eff, S2, MODE_15N);
					break;
				case R_13CR1:
					calc_R = SMF_R1(resid, &(resid->relaxation[i]), tau_eff, S2, MODE_13C);
					break;
				case R_13CR1p:
					calc_R = SMF_R2(resid, &(resid->relaxation[i]), tau_eff, S2, MODE_13C);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			chisq += ((pow(resid->relaxation[i].R - calc_R, 2)) / pow(resid->relaxation[i].Rerror, 2));
		}
		return chisq;
	} else if (model == MOD_EMF || model == MOD_EMFT || model == MOD_DEMF || model == MOD_DEMFT) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double S2s = opts[1];
		long double S2NHd;
		long double tauf = opts[2];
		long double taus_eff=0, tauf_eff=0, Eas=0, Eaf=0, S2f = resid->S2_dipolar / S2s;
		if (model == MOD_EMFT) {
			Eas = opts[3];
			Eaf = opts[4];
		} else if (model == MOD_DEMFT) {
			Eas = opts[4];
			Eaf = opts[5];
		} else {
			taus_eff = taus;
			tauf_eff = tauf;
		}
		if (model == MOD_DEMF || model == MOD_DEMFT) 
			S2f = opts[3];
		
		//if (!isnan(opts[2]))
		//	printf("\t\t%Le\n", opts[2]);
		if (taus < 0 || tauf < 0)
			chisq += 100000000;
		if (S2s < 0 || S2s > 1 || S2f < 0 || S2f > 1)
			chisq += 100000000;
		if (S2s < resid->S2_dipolar)
			chisq += 100000000;
		if (tauf > taus)
			chisq += 100000000;
		if (tauf > upper_lim_tf || taus > upper_lim_ts || tauf < lower_lim_tf)
			chisq += 100000000;
		int i;
		for (i = 0; i < resid->n_relaxation; i++) {
			if (resid->relaxation[i].R <= 0)
				continue;

			if (model == MOD_EMFT || model == MOD_DEMFT) {
				taus_eff = taus * expl(Eas / (RYD * resid->relaxation[i].T));
				tauf_eff = tauf * expl(Eaf / (RYD * resid->relaxation[i].T));
			}
			if (tauf_eff > taus_eff)
				chisq += 100000000;
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = EMF_R1(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, S2f, MODE_15N);
					break;
				case R_15NR1p:
					calc_R = EMF_R2(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, S2f, MODE_15N);
					break;
				case R_13CR1:
					calc_R = EMF_R1(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, S2f, MODE_13C);
					break;
				case R_13CR1p:
					calc_R = EMF_R2(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, S2f, MODE_13C);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			chisq += ((pow(resid->relaxation[i].R - calc_R, 2)) / pow(resid->relaxation[i].Rerror, 2));
		}
		if (model == MOD_DEMF || model == MOD_DEMFT ) { 
			S2NHd = S2s * S2f;
			chisq += ((pow(resid->S2_dipolar - S2NHd, 2)) / pow(0.01, 2));
		}
		return chisq;
	} else if (model == MOD_GAF || model == MOD_GAFT) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double tauf = opts[1];
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		long double sigf[3] = {opts[5], opts[6], opts[7]};
		long double taus_eff=0, tauf_eff=0, Eas=0, Eaf=0;
		if (model == MOD_GAFT) {
			Eas = opts[8];
			Eaf = opts[9];
		} else {
			taus_eff = taus;
			tauf_eff = tauf;
		}
		//if (!isnan(opts[2]))
		//	printf("\t\t%Le\n", opts[2]);
		int i;
		if (taus < 0 || tauf < 0)
			chisq += 100000000;
		for (i = 0; i < 3; i++) {
			if (sigs[i] < 0 || sigs[i] > 0.52360)
				chisq += 100000000;
			if (sigf[i] < 0 || sigf[i] > 0.52360)
				chisq += 100000000;
		}

		if (tauf > taus)
			chisq += 100000000;
		if (tauf > upper_lim_tf || taus > upper_lim_ts)
			chisq += 100000000;

		for (i = 0; i < resid->n_relaxation; i++) {
			if (resid->relaxation[i].R <= 0)
				continue;
			if (model == MOD_GAFT) {
				taus_eff = taus * expl(Eas / (RYD * resid->relaxation[i].T));
				tauf_eff = tauf * expl(Eaf / (RYD * resid->relaxation[i].T));
			}
			if (tauf_eff > taus_eff)
				chisq += 100000000;
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = GAF_15NR1(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				case R_15NR1p:
					calc_R = GAF_15NR2(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				case R_13CR1:
					calc_R = GAF_13CR1(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				case R_13CR1p:
					calc_R = GAF_13CR2(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			//printf("%f\n", calc_R);

			chisq += ((pow(resid->relaxation[i].R - calc_R, 2)) / pow(resid->relaxation[i].Rerror, 2));
		}

		double S2NHs, S2NHf, S2NHd;
		struct Orient *Os[] = {&(resid->orients[OR_NH])};
		double *S2s[] = {&S2NHs};
		double *S2f[] = {&S2NHf};
		GAF_S2(sigs, Os, Os, S2s, 1, MODE_REAL);
		GAF_S2(sigf, Os, Os, S2f, 1, MODE_REAL);
		
		S2NHd = S2NHs * S2NHf;
		chisq += ((pow(resid->S2_dipolar - S2NHd, 2)) / pow(0.01, 2));

		return chisq;
	} else {
		printf("Model not implemented yet\n");
		exit(-1);
	}
	return chisq; 
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
int back_calculate(long double * opts, struct Residue * resid, int model, char *filename) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */
	//if (opts[0] == -1)
	//	return -1;
	double calc_R = 0;
	FILE * fp;
	fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}


	if (model == MOD_SMF || model == MOD_SMFT) {
		/* Simple Model Free Analysis */
		long double tau = opts[0];
		long double S2 = opts[1];
		long double tau_eff=0, Ea=0;
		if (model == MOD_SMFT)
			Ea = opts[2];
		else
			tau_eff = tau;


		int i;
		for (i = 0; i < resid->n_relaxation; i++) {
			if (resid->relaxation[i].R <= 0)
                                continue;

			if (model == MOD_SMFT)
				tau_eff = tau * expl(Ea / (RYD * resid->relaxation[i].T));
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = SMF_R1(resid, &(resid->relaxation[i]), tau_eff, S2, MODE_15N);
					break;
				case R_15NR1p:
					calc_R = SMF_R2(resid, &(resid->relaxation[i]), tau_eff, S2, MODE_15N);
					break;
				case R_13CR1:
					calc_R = SMF_R1(resid, &(resid->relaxation[i]), tau_eff, S2, MODE_13C);
					break;
				case R_13CR1p:
					calc_R = SMF_R2(resid, &(resid->relaxation[i]), tau_eff, S2, MODE_13C);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			fprintf(fp, "%d\t%f\t%f\t%f", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);
			if (VERBOSE)
				fprintf(fp, "\t%f\t%f\t%f\t%f", resid->relaxation[i].field, resid->relaxation[i].wr, resid->relaxation[i].w1, resid->relaxation[i].T); 
			fprintf(fp, "\n");
		}
		fclose(fp);
		return 1;
	} else if (model == MOD_EMF || model == MOD_EMFT || model == MOD_DEMF || model == MOD_DEMFT) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double S2s = opts[1];
		long double tauf = opts[2];
		long double taus_eff=0, tauf_eff=0, Eas=0, Eaf=0, S2f = resid->S2_dipolar / S2s;
		if (model == MOD_EMFT) {
			Eas = opts[3];
			Eaf = opts[4];
		} else if (model == MOD_DEMFT) {
			Eas = opts[4];
			Eaf = opts[5];
		} else {
			taus_eff = taus;
			tauf_eff = tauf;
		}
		if (model == MOD_DEMF || model == MOD_DEMFT)
			S2f = opts[3];
		int i;
		//double EMF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf)
		for (i = 0; i < resid->n_relaxation; i++) {
			if (resid->relaxation[i].R <= 0)
                                continue;

			if (model == MOD_EMFT || model == MOD_DEMFT) {
				taus_eff = taus * expl(Eas / (RYD * resid->relaxation[i].T));
				tauf_eff = tauf * expl(Eaf / (RYD * resid->relaxation[i].T));
			}
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = EMF_R1(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, S2f, MODE_15N);
					break;
				case R_15NR1p:
					calc_R = EMF_R2(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, S2f, MODE_15N);
					break;
				case R_13CR1:
					calc_R = EMF_R1(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, S2f, MODE_13C);
					break;
				case R_13CR1p:
					calc_R = EMF_R2(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, S2f, MODE_13C);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			fprintf(fp, "%d\t%f\t%f\t%f", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);
			if (VERBOSE)
				fprintf(fp, "\t%f\t%f\t%f\t%f", resid->relaxation[i].field, resid->relaxation[i].wr, resid->relaxation[i].w1, resid->relaxation[i].T); 
			fprintf(fp, "\n");
		}
		fclose(fp);
		return 1;
	} else if (model == MOD_GAF || model == MOD_GAFT) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double tauf = opts[1];
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		long double sigf[3] = {opts[5], opts[6], opts[7]};
		long double taus_eff=0, tauf_eff=0, Eas=0, Eaf=0;
		if (model == MOD_GAFT) {
			Eas = opts[8];
			Eaf = opts[9];
		} else {
			taus_eff = taus;
			tauf_eff = tauf;
		}
		//if (!isnan(opts[2]))
		//	printf("\t\t%Le\n", opts[2]);
		int i;

		for (i = 0; i < resid->n_relaxation; i++) {
			if (resid->relaxation[i].R <= 0)
				continue;
			if (model == MOD_GAFT) {
				taus_eff = taus * expl(Eas / (RYD * resid->relaxation[i].T));
				tauf_eff = tauf * expl(Eaf / (RYD * resid->relaxation[i].T));
				//printf("%Le, %Le\n", taus_eff, tauf_eff);
			}

			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = GAF_15NR1(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				case R_15NR1p:
					calc_R = GAF_15NR2(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				case R_13CR1:
					calc_R = GAF_13CR1(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				case R_13CR1p:
					calc_R = GAF_13CR2(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			fprintf(fp, "%d\t%f\t%f\t%f", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);
			if (VERBOSE)
				fprintf(fp, "\t%f\t%f\t%f\t%f", resid->relaxation[i].field, resid->relaxation[i].wr, resid->relaxation[i].w1, resid->relaxation[i].T); 
			fprintf(fp, "\n");
		}
		fclose(fp);
		return 1;
	} else {
		printf("Model not implemented yet\n");
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
	return -1;
}
