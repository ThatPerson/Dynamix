#include <stdio.h>
#include <math.h>

double optimize_chisq(long double * opts, struct Residue * resid, int model) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */
	double calc_R = 0;
	double chisq = 0;
	long double upper_lim_tf = (long double) 1 * pow(10, -8);
	long double upper_lim_ts = (long double) 1 * pow(10, -5);
	//printf("%d %d: %d\n", model, MOD_SMFT, model == MOD_SMFT);
	if (model == MOD_SMF || model == MOD_SMFT) {
		/* Simple Model Free Analysis */
		long double tau = opts[0];
		long double S2 = opts[1];
		long double tau_eff;
		long double Ea;
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
	} else if (model == MOD_EMF || model == MOD_EMFT) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double S2s = opts[1];
		long double tauf = opts[2];
		long double taus_eff, tauf_eff, Eas, Eaf;
		if (model == MOD_EMFT) {
			Eas = opts[3];
			Eaf = opts[4];
		} else {
			taus_eff = taus;
			tauf_eff = tauf;
		}
		//if (!isnan(opts[2]))
		//	printf("\t\t%Le\n", opts[2]);
		if (taus < 0 || tauf < 0)
			chisq += 100000000;
		if (S2s < 0 || S2s > 1)
			chisq += 100000000;
		if (S2s < resid->S2_dipolar)
			chisq += 100000000;
		if (tauf > taus)
			chisq += 100000000;
		if (tauf > upper_lim_tf || taus > upper_lim_ts)
			chisq += 100000000;
		int i;
		for (i = 0; i < resid->n_relaxation; i++) {
			if (model == MOD_EMFT) {
				taus_eff = taus * expl(Eas / (RYD * resid->relaxation[i].T));
				tauf_eff = tauf * expl(Eaf / (RYD * resid->relaxation[i].T));
			}
			if (tauf_eff > taus_eff)
				chisq += 100000000;
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = EMF_R1(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, MODE_15N);
					break;
				case R_15NR1p:
					calc_R = EMF_R2(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, MODE_15N);
					break;
				case R_13CR1:
					calc_R = EMF_R1(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, MODE_13C);
					break;
				case R_13CR1p:
					calc_R = EMF_R2(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, MODE_13C);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			chisq += ((pow(resid->relaxation[i].R - calc_R, 2)) / pow(resid->relaxation[i].Rerror, 2));
		}
		return chisq;
	} else if (model == MOD_GAF || model == MOD_GAFT) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double tauf = opts[1];
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		long double sigf[3] = {opts[5], opts[6], opts[7]};
		long double taus_eff, tauf_eff, Eas, Eaf;
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
			if (sigs[i] < 0 || sigs[i] > 30)
				chisq += 100000000;
			if (sigf[i] < 0 || sigf[i] > 30)
				chisq += 100000000;
		}

		if (tauf > taus)
			chisq += 100000000;
		if (tauf > upper_lim_tf || taus > upper_lim_ts)
			chisq += 100000000;

		for (i = 0; i < resid->n_relaxation; i++) {
			if (model == MOD_EMFT) {
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
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			//printf("%f\n", calc_R);

			chisq += ((pow(resid->relaxation[i].R - calc_R, 2)) / pow(resid->relaxation[i].Rerror, 2));
		}
		return chisq;
	} else {
		printf("Model not implemented yet\n");
		exit(-1);
	}
	return -1;
}

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
		long double tau_eff, Ea;
		if (model == MOD_SMFT)
			Ea = opts[2];
		else
			tau_eff = tau;


		int i;
		for (i = 0; i < resid->n_relaxation; i++) {
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
			fprintf(fp, "%d\t%f\t%f\t%f\n", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);
		}
		fclose(fp);
		return 1;
	} else if (model == MOD_EMF || model == MOD_EMFT) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double S2s = opts[1];
		long double tauf = opts[2];
		long double taus_eff, tauf_eff, Eas, Eaf;
		if (model == MOD_EMFT) {
			Eas = opts[3];
			Eaf = opts[4];
		} else {
			taus_eff = taus;
			tauf_eff = tauf;
		}



		int i;
		//double EMF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf)
		for (i = 0; i < resid->n_relaxation; i++) {
			if (model == MOD_EMFT) {
				taus_eff = taus * expl(Eas / (RYD * resid->relaxation[i].T));
				tauf_eff = tauf * expl(Eaf / (RYD * resid->relaxation[i].T));
			}
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = EMF_R1(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, MODE_15N);
					break;
				case R_15NR1p:
					calc_R = EMF_R2(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, MODE_15N);
					break;
				case R_13CR1:
					calc_R = EMF_R1(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, MODE_13C);
					break;
				case R_13CR1p:
					calc_R = EMF_R2(resid, &(resid->relaxation[i]), taus_eff, S2s, tauf_eff, MODE_13C);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			fprintf(fp, "%d\t%f\t%f\t%f\n", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);
		}
		fclose(fp);
		return 1;
	} else if (model == MOD_GAF || model == MOD_GAFT) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double tauf = opts[1];
		long double sigs[3] = {opts[2], opts[3], opts[4]};
		long double sigf[3] = {opts[5], opts[6], opts[7]};
		long double taus_eff, tauf_eff, Eas, Eaf;
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
			if (model == MOD_EMFT) {
				taus_eff = taus * expl(Eas / (RYD * resid->relaxation[i].T));
				tauf_eff = tauf * expl(Eaf / (RYD * resid->relaxation[i].T));
			}

			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = GAF_15NR1(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				case R_15NR1p:
					calc_R = GAF_15NR2(resid, &(resid->relaxation[i]), taus_eff, tauf_eff, sigs, sigf);
					break;
				default:
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			fprintf(fp, "%d\t%f\t%f\t%f\n", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);
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
