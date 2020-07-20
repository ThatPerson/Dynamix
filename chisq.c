#include <stdio.h>
#include <math.h>

double optimize_chisq(long double * opts, struct Residue * resid, int model) {
	/* opts is a pointer to an array containing;
	 * 
	 *  for SMF, [tau, S2]
	 */
	double calc_R = 0;
	double chisq = 0;
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
	} else {
		printf("Model not implemented yet\n");
		exit(-1);
	}
	return -1;
}
