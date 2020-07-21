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
	if (model == MOD_SMF) {
		/* Simple Model Free Analysis */
		long double tau = opts[0];
		long double S2 = opts[1];
		if (tau < 0)
			chisq += 100000000;
		if (S2 < 0 || S2 > 1)
			chisq += 100000000;
		int i;
		for (i = 0; i < resid->n_relaxation; i++) {
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = SMF_15NR1(resid, &(resid->relaxation[i]), tau, S2);
					break;
				case R_15NR1p:
					calc_R = SMF_15NR2(resid, &(resid->relaxation[i]), tau, S2);
					break;
				default: 
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			chisq += ((pow(resid->relaxation[i].R - calc_R, 2)) / pow(resid->relaxation[i].Rerror, 2));
		}
		return chisq;
	} else if (model == MOD_EMF) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double S2s = opts[1];
		long double tauf = opts[2];
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
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = EMF_15NR1(resid, &(resid->relaxation[i]), taus, S2s, tauf);
					break;
				case R_15NR1p:
					calc_R = EMF_15NR2(resid, &(resid->relaxation[i]), taus, S2s, tauf);
					break;
				default: 
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
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
	double calc_R = 0;
	FILE * fp;
	fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	
	
	if (model == MOD_SMF) {
		/* Simple Model Free Analysis */
		long double tau = opts[0];
		long double S2 = opts[1];

		int i;
		for (i = 0; i < resid->n_relaxation; i++) {
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = SMF_15NR1(resid, &(resid->relaxation[i]), tau, S2);
					break;
				case R_15NR1p:
					calc_R = SMF_15NR2(resid, &(resid->relaxation[i]), tau, S2);
					break;
				default: 
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			fprintf(fp, "%d\t%f\t%f\t%f\n", i, calc_R, resid->relaxation[i].R, resid->relaxation[i].Rerror); 
		}
		fclose(fp);
		return 1;
	} else if (model == MOD_EMF) {
		/* Extended Model Free Analysis */
		long double taus = opts[0];
		long double S2s = opts[1];
		long double tauf = opts[2];
		int i;
		//double EMF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf) 
		for (i = 0; i < resid->n_relaxation; i++) {
			switch (resid->relaxation[i].type) {
				case R_15NR1:
					calc_R = EMF_15NR1(resid, &(resid->relaxation[i]), taus, S2s, tauf);
					break;
				case R_15NR1p:
					calc_R = EMF_15NR2(resid, &(resid->relaxation[i]), taus, S2s, tauf);
					break;
				default: 
					printf("Unknown relaxation type: %d\n", resid->relaxation[i].type);
					exit(0);
					break;
			}
			fprintf(fp, "%d\t%f\t%f\t%f\n", i, calc_R, resid->relaxation[i].R, resid->relaxation[i].Rerror); 
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
