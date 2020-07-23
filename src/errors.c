#include <stdio.h>
#include <math.h>

long double uniform_rand(void) {
	return ((long double) rand() + 1.) / ((long double) RAND_MAX + 1.);
}

float norm_rand(float mean, float std) {
	// Box-Muller method
	float rnd1, rnd2;
	rnd1 = (float) uniform_rand();
	rnd2 = (float) uniform_rand();
	float unadj = sqrt(-2 * log(rnd1)) * cos(2 * M_PI * rnd2);
	
	//printf("%f, %f, %f\n", rnd1, rnd2, mean + std*unadj);
	return mean + std * unadj;
}

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


void * calc_errors(void *input) {
	int i = ((struct rrargs*)input)->i;
	printf("\tThread %d alive...\n", i+1);
	struct Residue * resid = ((struct rrargs*)input)->resid;
	int model = ((struct rrargs*)input)->model;
	int n_iter = ((struct rrargs*)input)->n_iter;
	//char outputdir[255];
	//strcpy(outputdir, ((struct rrargs*)input)->outputdir);


	int l, k, params = 0;
	switch (model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT: params= 5; break;
		case MOD_SMFT: params= 3; break;
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
	double val;
	for (l = 0; l < n_iter; l++) {
		if (resid->ignore == 1) {
			return NULL;
		}
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
		
		val = simplex(optimize_chisq, opts, params, 1.0e-16, 1, resid, model);
		
		// The actual value is more or less irrelevant, we're just interested in the values now in 'opts'
		//resid->error_params[l] = (long double *) malloc (sizeof(long double) * params);
		for (k = 0; k < params; k++) {
			resid->error_params[k][l] = opts[k];
		}
		
		
		/* Return our pointers to how they were before... */
		free(resid->relaxation);
		resid->relaxation = NULL;
		resid->relaxation = resid->temp_relaxation;
		resid->temp_relaxation = NULL;
	}
	free(opts);
	return NULL;
}
