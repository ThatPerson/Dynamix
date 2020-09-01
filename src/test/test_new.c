/**
 * @file main.c
 */

#include <stdio.h>
#include "datatypes.c"
#include "read_data.c"
#include "models.c"
#include "chisq.c"
#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
#include "errors.c"
#include <time.h>
#include <pthread.h>

/**
 * Operates residue optimization. Generates random parameter guesses and passes these to the simplex function.
 * @param input
 *  Pointer to rrarg containing thread information
 */
void * run_residue(void *input) {
	int i = ((struct rrargs*)input)->i;
	printf("\tThread %d alive...\n", i + 1);
	struct Residue * resid = ((struct rrargs*)input)->resid;
	int model = ((struct rrargs*)input)->model;
	int n_iter = ((struct rrargs*)input)->n_iter;
	char outputdir[255];
	strcpy(outputdir, ((struct rrargs*)input)->outputdir);
	FILE *fp;
	char filename[300];
	sprintf(filename, "%s/residue_%d.dat", outputdir, i+1);
	fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return NULL;
	}

	//printf("RESIDUE %d\n", i+1);
	//printf("Number of relaxations: %d\n", resid->n_relaxation);
	int l, k, params = 0;
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
	resid->parameters = (long double *) malloc (sizeof(long double) * params);
	resid->min_val = MIN_VAL;

	for (l = 0; l < n_iter; l++) {
		if (resid->ignore == 1) {
			fprintf(fp, "%d, %f, -1, -1\n", i+1, 1000.);
			fclose(fp);
			free(opts);
			return NULL;
		}
		//long double opts[20] = {0, 0};
		
		/**
		 * SMF parameters are \n 
		 *   [0] tau\n 
		 *   [1] S2\n 
		 * SMFT parameters;\n 
		 *   [0] tau\n 
		 *   [1] S2\n 
		 *   [2] Ea\n 
		 * EMF parameters\n 
		 *   [0] tau slow\n 
		 *   [1] S2 slow\n 
		 *   [2] tau fast\n 
		 *   NOTE: The fast order parameter is calculated as S2_dipolar/S2s\n 
		 * EMFT parameters\n 
		 *   [0] tau slow\n 
		 *   [1] S2 slow\n 
		 *   [2] tau fast\n 
		 *   [3] activation energy for slow motion\n 
		 *   [4] activation energy for fast motion\n 
		 * DEMF parameters\n
		 *   [0] tau slow\n
		 *   [1] S2 slow\n
		 *   [2] tau fast\n
		 *   [3] S2 fast\n
		 * DEMFT parameters\n
		 *   [0] tau slow\n
		 *   [1] S2 slow\n
		 *   [2] tau fast\n
		 *   [3] S2 fast\n
		 *   [4] activation energy for slow motion\n
		 *   [5] activation energy for fast motion\n
		 * GAF parameters\n 
		 *   [0] tau slow\n 
		 *   [1] tau fast\n 
		 *   [2-4] alpha, beta, gamma deflections for slow motions\n 
		 *   [5-7] alpha, beta, gamma deflections for fast motions\n 
		 * GAFT parameters\n 
		 *   [0] tau slow\n 
		 *   [1] tau fast\n 
		 *   [2-4] alpha, beta, gamma deflections for slow motions\n 
		 *   [5-7] alpha, beta, gamma deflections for fast motions\n 
		 *   [8] activation energy for slow motion\n 
		 *   [9] activation energy for fast motion\n 
		 */
		if (model == MOD_SMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
		} else if (model == MOD_EMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
			//printf("RUN: %f, %Le, %Le, %Le\n", resid->S2_dipolar, opts[0], opts[1], opts[2]);
		} else if (model == MOD_EMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
			opts[3] = (rand()%60000)/1.;
			opts[4] = (rand()%60000)/1.;
			//printf("RUN: %f, %Le, %Le, %Le\n", resid->S2_dipolar, opts[0], opts[1], opts[2]);
		} else if (model == MOD_DEMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
			opts[3] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.);
		} else if (model == MOD_DEMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
			opts[3] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.);
			opts[4] = (rand()%60000)/1.;
			opts[5] = (rand()%60000)/1.;
		} else if (model == MOD_SMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
			opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
			opts[2] = (rand()%60000)/1.;
		} else if (model == MOD_GAF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
			for (k = 2; k <= 7; k++) {
				// 15 degrees = 0.26180 radians
				opts[k] = ((rand () % 250)/1000.);
				//printf("%d %Lf\n", k, opts[k] * (180. / M_PI));
			}
		} else if (model == MOD_GAFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
			opts[1] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
			for (k = 2; k <= 7; k++) {
				// 15 degrees = 0.26180 radians
				opts[k] = ((rand () % 250)/1000.);
			}
			opts[8] = (rand()%60000)/1.;
			opts[9] = (rand()%60000)/1.;
		}

		//printf("%Le, %Le\n", opts[0], opts[1]);
		double val = simplex(optimize_chisq, opts, params, 1.0e-16, 1, resid, model);
		if (val >= 1000000) {
			val = -1;
			for (k = 0; k < params; k++) {
				opts[k] = -1;
			}
		}
		fprintf(fp, "%d\t%f", i+1, val);
		for (k = 0; k < params; k++) {
			fprintf(fp, "\t%Le", opts[k]);
		}

		if (val < resid->min_val && val != -1) {
			//printf("New lowest %f\n", val);
			resid->min_val = val;
			for (k = 0; k < params; k++) {
				resid->parameters[k] = opts[k];
			}
		}
		fprintf(fp, "\n");
	}
	free(opts);
	fclose(fp);
	return NULL;
}



/**
 * Start function. Initialises parameters, loads files, spawns threads and outputs data.
 * @opts filename
 *  Takes input as path to file containing System definition (*dx file)
 * @opts -e
 *  Enables error mode
 */
int main(int argc, char * argv[]) {

	/* Initialisation */
	srand(time(NULL));
	initialise_dwig();

	char system_file[255] = "";
	int i;
	int err_mod = 0;
	for (i = 1; i < argc; i++) {
		//printf("%s\n", argv[i]);
		if (strcmp(argv[i], "-e") == 0)
			err_mod = 1;
		else
			strcpy(system_file, argv[i]);
	}

	if (strcmp(system_file, "") == 0){
		printf("Please provide system file.\n");
		exit(-1);
	}



	struct Model m;
	int ret = read_system_file(system_file, &m);
	
	
	
	m.error_mode = err_mod;
	if (m.error_mode == 1 && m.n_error_iter < 0) {
		printf("Please provide number of error iterations\n");
		ret = -1;
	}

	if (ret == -1) {
		printf("Error found, crashing peacefully...\n");
		exit(-1);
	}

	printf("%d\n", ret);
	char filename[300];
	sprintf(filename, "%s/model.txt", m.outputdir);
	print_system(&m, filename);


	//double optimize_chisq(long double * opts, struct Residue * resid, int model) {

	/* Test Case 1. Residue 14 appears to give 20677.643567 chisq for _any_ params. 
	  - this is the NaN error case (eg it's just the sum of values over error squares)*/
	//long double opts_o[10] = {4.000000e-08, 3.100000e-12, 5.900000e-02, 1.520000e-01, 1.970000e-01, 1.370000e-01, 1.370000e-01, 1.990000e-01, 3.956100e+04, 1.831000e+04};
	//long double opts_o[5] = {1.372975e+00, 7.989054e-01, 5.587675e-04, 3.777893e+04, 2.549346e+04};
	//long double opts_t[10] = {5.800000e-07, 3.500000e-12, 8.500000e-02, 1.760000e-01, 8.700000e-02, 8.000000e-03, 2.140000e-01, 2.180000e-01, 5.341000e+04, 4.957900e+04};
	//long double opts_o[3] = {7.558008e+00, 1.000000e+00, 4.937203e-02};
	//long double opts_t[3] = {8.436836e+00, 8.529099e-01, 1.000000e-04};
	long double opts_o[] = {4.900000e+01, 9.600000e-2, 1.100000e-02, 1.160000e-01, 7.600000e-02, 2.340000e-01, 1.990000e-01, 2.190000e-01, 4.936500e+04, 4.540900e+04};
	long double opts_t[] = {2.300000e+01, 8.000000e-3, 8.300000e-02, 2.320000e-01, 9.100000e-02, 8.400000e-02, 1.770000e-01, 1.990000e-01, 3.513400e+04, 9.790000e+03};
	double result_o = optimize_chisq(opts_o, &(m.residues[35]), MOD_GAFT);
	double result_t = 0;
	result_t = optimize_chisq(opts_t, &(m.residues[35]), MOD_GAFT);
	
	// 1181.839707     8.436836e+00    8.529099e-01    1.000000e-04
	printf("%f %f\n", result_o, result_t);
	
	
	
	/*fclose(fp);
	if (m.error_mode == 1)
		fclose(ep);
	if (m.model == MOD_GAF || m.model == MOD_GAFT)
		fclose(gaf);*/
	free_all(&m);
	return 1;
}
