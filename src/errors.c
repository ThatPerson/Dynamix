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
	*std = sqrtl(m2 - powl(*mean, 2));

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
	unsigned int n_iter = ((struct rrargs*)input)->n_iter;
	//double optim = resid->min_val;
	//char outputdir[255];
	//strcpy(outputdir, ((struct rrargs*)input)->outputdir);


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
		default: params = 0; break;
	}
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

		for (k = 0; k < resid->n_relaxation; k++) {
			resid->relaxation[k].field = resid->temp_relaxation[k].field;
			resid->relaxation[k].wr = resid->temp_relaxation[k].wr;
			resid->relaxation[k].w1 = resid->temp_relaxation[k].w1;
			resid->relaxation[k].type = resid->temp_relaxation[k].type;
			resid->relaxation[k].T = resid->temp_relaxation[k].T;
			resid->relaxation[k].R = norm_rand(resid->temp_relaxation[k].R, (resid->temp_relaxation[k].Rerror)/2.);
			//printf("%f, %f -> %f \n", resid->temp_relaxation[k].R, resid->temp_relaxation[k].Rerror, resid->relaxation[k].R);
			/* Rerror is Monte-Carlo calculated 2 standard deviations; norm_rand takes 1 std */
			resid->relaxation[k].Rerror = -1.;
		}
		for (k = 0; k < params; k++) {
			opts[k] = resid->parameters[k];
		}
		simplex(optimize_chisq, opts, params, 1.0e-16, 1, resid, model);

		// The actual value is more or less irrelevant, we're just interested in the values now in 'opts'
		//resid->error_params[l] = (long double *) malloc (sizeof(long double) * params);
		for (k = 0; k < params; k++) {
			resid->error_params[k][p] = opts[k];
		}
		p++;


		/* Return our pointers to how they were before... */
		free(resid->relaxation);
		resid->relaxation = NULL;
		resid->relaxation = resid->temp_relaxation;
		resid->temp_relaxation = NULL;
	}
	resid->error_calcs = p;
	free(opts);
	return NULL;
}
