/**
 * @file errors.c
 */

#include <stdio.h>
#include <math.h>

/** 
 * Generates uniform random long double from 0 to 1 inclusive.
 * @return long double
 *  Long double containing uniform random number.
 */
long double uniform_rand(void) {
	return ((long double) rand() + 1.) / ((long double) RAND_MAX + 1.);
}

/**
 * Uses Box-Muller method to generate a normally distributed random number.
 * @param mean
 *  Mean of normal distribution to select random variable from
 * @param std
 *  Standard deviation of normal distribution from which random variable selected
 * @return float
 *  Returns normally distributed random number
 */
double norm_rand(double mean, double std) {
	// Box-Muller method
	double rnd1, rnd2;
	rnd1 = (double) uniform_rand();
	rnd2 = (double) uniform_rand();
	double unadj = sqrt(-2 * log(rnd1)) * cos(2 * M_PI * rnd2);

	//printf("%f, %f, %f\n", rnd1, rnd2, mean + std*unadj);
	return mean + std * unadj;
}

/**
 * Calculates mean and standard deviation of values contained in array, then puts these into the given pointers.
 * @param vals
 *  Pointer to array of values to take statistics of
 * @param length
 *  Length of vals
 * @param mean
 *  Pointer to long double to contain mean
 * @param std
 *  Pointer to long double to contain standard deviation
 */
void calc_statistics(long double * vals, int length, long double * mean, long double * std) {
	int k;
	*mean = 0;
	long double m2 = 0;
	for (k = 0; k < length; k++) {
		*mean += vals[k];
		m2 += powl(vals[k], 2);
	}
	*mean = *mean / ((long double) length);
	m2 = m2 / ((long double) length);
	long double mdiff = m2 - powl(*mean, 2.);
	if (mdiff < 0 && fabsl(mdiff) < 0.000001) {
		ERROR("mdiff = %Lf, making negative.", mdiff);
		// then all have converged to the same point and the std is 0. 
		// but because floating points are ew, it will give 'nan' for sqrtl(-0.0000)
		// so we just invert the sign. It's not pretty but ah well.
		mdiff = -mdiff;
	}
	*std = sqrtl(mdiff);

}

/**
 * Calculates errors for a given residue.
 * Error calculation is done by varying the relaxation values according to a normal distribution with mean (R) and standard deviation (Rerror/2).
 * Then simplex optimization is performed, and the newly optimized parameters stored. Then statistics of these back optimized parameters are taken,
 * with the errors being the standard deviation of these.
 * @param input
 *  rrarg struct containing the residue under consideration.
 */
void * calc_errors(void *input) {
	unsigned int i = ((struct rrargs*)input)->i;
	printf("\tThread %d alive...\n", i+1);
	struct Residue * resid = ((struct rrargs*)input)->resid;
	unsigned int model = ((struct rrargs*)input)->model;
	unsigned int or_variation = ((struct rrargs*)input)->or_variation;
	unsigned int n_iter = ((struct rrargs*)input)->n_iter;
	char outputdir[255];
	strcpy(outputdir, ((struct rrargs*)input)->outputdir);

	//double optim = resid->min_val;
	//char outputdir[255];
	//strcpy(outputdir, ((struct rrargs*)input)->outputdir);
	FILE *errp;
	char filename[300];
	sprintf(filename, "%s/errors_%d.dat", outputdir, i+1);
	errp = fopen(filename, "w");
	if (errp == NULL) {
		printf("%s not found.\n", filename);
		return NULL;
	}

	unsigned int l, k, params = 0;
	switch (model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT: params= 5; break;
		case MOD_SMFT: params= 3; break;
		case MOD_DEMF: params = 4; break;
		case MOD_DEMFT: params= 6; break;
		case MOD_GAF: params = 8; break;
		case MOD_GAFT:params = 10; break;
		case MOD_EGAF: params = 6; break;
		case MOD_EGAFT: params = 8; break;
		default: params = 0; break;
	}
	if (or_variation == VARIANT_A)
                params += 3; // alpha, beta, gamma

	//printf("%d\n", params);
	long double *opts;
	opts = (long double *) malloc (sizeof(long double) * params);
	resid->errors_mean = (long double *) malloc (sizeof(long double) * params);
	resid->errors_std = (long double *) malloc(sizeof(long double) * params);
	resid->error_params = (long double **) malloc (sizeof(long double *) * params);
	for (k = 0; k < params; k++) {
		resid->error_params[k] = (long double *) malloc (sizeof(long double) * n_iter);
	}
	//resid->min_val = MIN_VAL;*/
	//double val = 0;
	int p = 0;
	int ignore = -1;
	double temp_R = 0;

	// make backups of order parameters.
	resid->S2NHb = resid->S2NH;
	resid->S2CHb = resid->S2CH;
	resid->S2CNb = resid->S2CN;

	for (l = 0; l < n_iter; l++) {
		if (resid->ignore == 1) {
			return NULL;
		}
		
		/* As the simplex optimization looks in the relaxation pointer, we temporarily
		 * store the relaxation pointer elsewhere while we perform the optimization
		 */
		resid->temp_relaxation = (resid->relaxation);
		resid->relaxation = NULL;
		resid->relaxation = (struct Relaxation *) malloc(sizeof(struct Relaxation) * resid->lim_relaxation);

		//if (i == 2)
		//	printf("Residue %d iteration %d: %Lf\n ", i, l, opts[0]);

		for (k = 0; k < resid->n_relaxation; k++) {
			resid->relaxation[k].field = resid->temp_relaxation[k].field;
			resid->relaxation[k].wr = resid->temp_relaxation[k].wr;
			resid->relaxation[k].w1 = resid->temp_relaxation[k].w1;
			resid->relaxation[k].type = resid->temp_relaxation[k].type;
			resid->relaxation[k].T = resid->temp_relaxation[k].T;

			temp_R = back_calc(resid->parameters, resid, &(resid->temp_relaxation[k]), model, &ignore);
			//if (i == 2 && k == 2)
			//	printf("Iteration %d, R = %f\n", l, temp_R); // should be the same each time
			resid->relaxation[k].R = norm_rand(temp_R, (resid->temp_relaxation[k].Rerror/2.));


			LOG("%d %d %d prior: %f, backcalc: %f, new: %f, error: %f", i, l, k, resid->temp_relaxation[k].R, temp_R, resid->relaxation[k].R, resid->temp_relaxation[k].Rerror/2.);

			//resid->relaxation[k].R = norm_rand(resid->temp_relaxation[k].R, (resid->temp_relaxation[k].Rerror)/2.);
			//printf("%f, %f -> %f \n", resid->temp_relaxation[k].R, resid->temp_relaxation[k].Rerror, resid->relaxation[k].R);
			/* Rerror is Monte-Carlo calculated 2 standard deviations; norm_rand takes 1 std */
			resid->relaxation[k].Rerror = resid->temp_relaxation[k].Rerror;
		}

		resid->S2NH = norm_rand(resid->S2NHb, (resid->S2NHe / 2.)); // divided by two to get 1 standard deviation
		resid->S2CH = norm_rand(resid->S2CHb, (resid->S2CHe / 2.));
		resid->S2CN = norm_rand(resid->S2CNb, (resid->S2CNe / 2.));

		for (k = 0; k < params; k++) {
			opts[k] = resid->parameters[k];
		}
		double min = simplex(optimize_chisq, opts, params, 1.0e-16, 1, resid, model, or_variation);
		LOG("%d %d min = %f; orig = %f", i, l, min, resid->min_val);
		// The actual value is more or less irrelevant, we're just interested in the values now in 'opts'
		//resid->error_params[l] = (long double *) malloc (sizeof(long double) * params);

		fprintf(errp, "%d\t%f", l+1, min);
		for (k = 0; k < params; k++) {
			resid->error_params[k][p] = opts[k];
			fprintf(errp, "\t%Le", opts[k]);
		}
		fprintf(errp, "\n");
		p++;


		/* Return our pointers to how they were before... */
		free(resid->relaxation);
		resid->relaxation = NULL;
		resid->relaxation = resid->temp_relaxation;
		resid->temp_relaxation = NULL;
	}

	// restore of order parameters.
	resid->S2NH = resid->S2NHb;
	resid->S2CH = resid->S2CHb;
	resid->S2CN = resid->S2CNb;

	resid->error_calcs = p;
	fclose(errp);
	free(opts);
	return NULL;
}
