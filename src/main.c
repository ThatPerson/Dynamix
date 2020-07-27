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
	char line[255];
	size_t len=255;
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
			opts[0] = ((rand() % 100)/100.) * powl(10, -8);
			opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
		} else if (model == MOD_EMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8);
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -11);
			//printf("RUN: %f, %Le, %Le, %Le\n", resid->S2_dipolar, opts[0], opts[1], opts[2]);
		} else if (model == MOD_EMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15);
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -20);
			opts[3] = (rand()%60000)/1.;
			opts[4] = (rand()%60000)/1.;
			//printf("RUN: %f, %Le, %Le, %Le\n", resid->S2_dipolar, opts[0], opts[1], opts[2]);
		} else if (model == MOD_DEMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8);
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -11);
			opts[3] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.);
		} else if (model == MOD_DEMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15);
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -20);
			opts[3] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.);
			opts[4] = (rand()%60000)/1.;
			opts[5] = (rand()%60000)/1.;
		} else if (model == MOD_SMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15);
			opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
			opts[2] = (rand()%60000)/1.;
		} else if (model == MOD_GAF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8);
			opts[1] = ((rand() % 100)/100.) * powl(10, -11);
			for (k = 2; k <= 7; k++) {
				// 15 degrees = 0.26180 radians
				opts[k] = ((rand () % 250)/1000.);
				//printf("%d %Lf\n", k, opts[k] * (180. / M_PI));
			}
		} else if (model == MOD_GAFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15);
			opts[1] = ((rand() % 100)/100.) * powl(10, -20);
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
	//print_system(&m, filename);

	/*int AS;
	long double sA, sB, sG;


	for (AS = 0; AS < 14; AS++) {
		sA = 2.591225;
		sB = 1.511158;
		sG = 0.684270;


		//printf("\n\n\nPhi %f, Theta %f\n", m.residues[28].orients[AS].phi, m.residues[28].orients[AS].theta);
		//printf("\nGAF_ord_paramTT(%Lf, %Lf, %Lf, %f, %f, %f, %f)\n\n", sA, sB, sG, m.residues[28].orients[AS].phi, m.residues[28].orients[AS].theta, m.residues[28].orients[AS].phi, m.residues[28].orients[AS].theta);

		//printf("\tsA: %Lf\n\tsB: %Lf\n\tsG: %Lf\n", sA, sB, sG);
		printf("\nGAF_ord_paramTT(%Lf, %Lf, %Lf, %f, %f, %f, %f)\n\n", sA, sB, sG, m.residues[28].orients[AS].phi, m.residues[28].orients[AS].theta, m.residues[28].orients[AS].phi, m.residues[28].orients[AS].theta);
		printf("\t\t... S2 = %0.36Lf\n", GAF_S2(sA, sB, sG, &(m.residues[28].orients[AS]), &(m.residues[28].orients[AS]), MODE_REAL));

		printf("\n\n\n");

	}*/


	pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * m.nthreads);
	pthread_attr_t threadattr;
	pthread_attr_init(&threadattr);
	pthread_attr_setstacksize(&threadattr, THREAD_STACK);
	int rc;

	int current_residue = 0;
	int n_spawns = 0;
	/* if we have 56 residues and 4 threads then we need
	 * 56 / 4 spawn events (= 14). Add 1 in case (eg for 57).
	 * Then loop over threads, increment current_residue and assign pointers.
	 */
	n_spawns = m.n_residues / m.nthreads;
	n_spawns ++;
	//printf("%d spawns\n", n_spawns);

	struct rrargs * RRA = (struct rrargs *)malloc(sizeof(struct rrargs) * m.nthreads);
	printf("Optimizing...\n");
	/* spawn the threads */
	int l = 0;
	for (l = 0; l < n_spawns; l++) {
		for (i=0; i<m.nthreads; ++i) {
			RRA[i].i = current_residue + i;
			if (current_residue + i >= m.n_residues)
				continue;
			RRA[i].resid = &(m.residues[current_residue + i]);
			RRA[i].model = m.model;
			RRA[i].n_iter = m.n_iter;
			strcpy(RRA[i].outputdir, m.outputdir);
			//printf("spawning thread %d (residue %d)\n", i, current_residue + i);
			rc = pthread_create(&threads[i], &threadattr, run_residue, (void *) &RRA[i]);
			if (rc != 0) {
				printf("Failed to spawn thread %d. Crashing gracefully...\n", current_residue + i);
				free_all(&m);
				free(RRA);
				free(threads);
				exit(-1);
			}
		}
		current_residue += m.nthreads;


		for (i=0; i<m.nthreads; ++i) {
			rc = pthread_join(threads[i], NULL);
		}
	}

	FILE * fp;
	sprintf(filename, "%s/final.dat", m.outputdir);
	fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		free_all(&m);
		return -1;
	}
	int params;
	FILE * ep = NULL;

	switch (m.model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT:params = 5; break;
		case MOD_SMFT:params = 3; break;
		case MOD_DEMF: params = 4; break;
		case MOD_DEMFT: params= 6; break;
		case MOD_GAF: params = 8; break;
		case MOD_GAFT:params = 10; break;
		default: params = 0; break;
	}

	if (m.error_mode == 1) {
		sprintf(filename, "%s/errors.dat", m.outputdir);
		ep = fopen(filename, "w");
		if (ep == NULL) {
			printf("%s not found.\n", filename);
			free_all(&m);
			return -1;
		}



		printf("Calculating Errors...\n");
		current_residue = 0;
		for (l = 0; l < n_spawns; l++) {
			for (i=0; i<m.nthreads; ++i) {
				RRA[i].i = current_residue + i;
				if (current_residue + i >= m.n_residues)
					continue;
				RRA[i].resid = &(m.residues[current_residue + i]);
				RRA[i].model = m.model;
				RRA[i].n_iter = m.n_error_iter;
				//strcpy(RRA[i].outputdir, m.outputdir);
				//printf("spawning thread %d (residue %d)\n", i, current_residue + i);
				rc = pthread_create(&threads[i], &threadattr, calc_errors, (void *) &RRA[i]);
				if (rc != 0) {
					printf("Failed to spawn thread %d. Crashing gracefully...\n", current_residue + i);
					free_all(&m);
					exit(-1);
				}
			}
			current_residue += m.nthreads;


			for (i=0; i<m.nthreads; ++i) {
				rc = pthread_join(threads[i], NULL);
			}
		}
	}

	free(RRA);
	free(threads);

	int k;
	FILE * gaf;
	if (m.model == MOD_GAF || m.model == MOD_GAFT) {
		sprintf(filename, "%s/gaf.dat", m.outputdir);
		gaf = fopen(filename, "w");
		if (gaf == NULL) {
			printf("%s not found.\n", filename);
			free_all(&m);
			return -1;
		}
	}

	printf("Outputting Files...\n");
	long double c = -1;
	char file[300];
	for (l = 0; l < m.n_residues; l++) {
		if (m.error_mode == 1) {
			for (k = 0; k < params; k++) {
				calc_statistics(m.residues[l].error_params[k], m.residues[l].error_calcs, &(m.residues[l].errors_mean[k]), &(m.residues[l].errors_std[k]));
			}
		}
		if (m.residues[l].min_val == MIN_VAL) {
			m.residues[l].min_val = -1.;
			for (i = 0; i < params; i++) {
				m.residues[l].parameters[i] = -1.;
			}
		}
		fprintf(fp, "%d\t%f\t%f", l+1, m.residues[l].S2_dipolar, m.residues[l].min_val);
		if (m.error_mode == 1) {
			fprintf(ep, "%d, %f, %f", l+1, m.residues[l].S2_dipolar, m.residues[l].min_val);
		}
		for (i = 0; i < params; i++) {
			fprintf(fp, "\t%Le", m.residues[l].parameters[i]);
			if (m.error_mode == 1) {
				fprintf(ep, ", %Le, %Le", m.residues[l].parameters[i], 2 * m.residues[l].errors_std[i]);
			}
			/* WARNING: I'm printing the actual minimized parameters with the errors from calculation.
			 * The error_means are generally not minimal. */
		}

		fprintf(fp, "\t%d\n", m.residues[l].error_calcs);
		if (m.error_mode == 1) {
			fprintf(ep, "\n");
		}
		sprintf(file, "%s/backcalc_%d.dat", m.outputdir, l+1);
		back_calculate((m.residues[l].parameters), &(m.residues[l]), m.model, file);
		
		if (m.model == MOD_GAF || m.model == MOD_GAFT) {
			// Print out 'effective S2' values.
			double S2slow, S2fast;
			double *S2[] = {&S2slow};
			/* Approximate as just the S2NHs and S2NHf */
			struct Orient *As[] = {&(m.residues[l].orients[OR_NH])};
			long double sigs[3] = {m.residues[l].parameters[2], m.residues[l].parameters[3], m.residues[l].parameters[4]};
			long double sigf[3] = {m.residues[l].parameters[5], m.residues[l].parameters[6], m.residues[l].parameters[7]};
			GAF_S2(sigs, As, As, S2, 1, MODE_REAL);
			S2[0] = &S2fast;
			GAF_S2(sigf, As, As, S2, 1, MODE_REAL);
			fprintf(gaf, "%d, %Le, %f, %Le, %f\n", l+1, m.residues[l].parameters[0], S2slow, m.residues[l].parameters[1], S2fast);
		}
	}


	fclose(fp);
	if (m.error_mode == 1)
		fclose(ep);
	if (m.model == MOD_GAF || m.model == MOD_GAFT)
		fclose(gaf);
	free_all(&m);
	return 1;
}
