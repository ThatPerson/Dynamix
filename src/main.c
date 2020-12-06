/**
 * @file main.c
 */

#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include "datatypes.h"
#include "read_data.h"

#include "models/smf.h"
#include "models/emf.h"
#include "models/gaf.h"
#include "models/egaf.h"
#include "chisq.h"
#include "crosen.h" // implementation of Nelder-Mead simplex algorithm
#include "errors.h"





#include "verification.h"


void * run_residue(void *input);

/**
 * Operates residue optimization. Generates random parameter guesses and passes these to the simplex function.
 * @param input
 *  Pointer to rrarg containing thread information
 */
void * run_residue(void *input) {
	unsigned int i = ((struct rrargs*)input)->i;
	printf("\tThread %d alive...\n", i + 1);
	struct Residue * resid = ((struct rrargs*)input)->resid;
	struct Model * m = ((struct rrargs*)input)->model;

	unsigned int model = m->model;
	unsigned int n_iter = ((struct rrargs*)input)->n_iter;
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
	unsigned int l, k;
	unsigned int params = 0;
	params = m->params;
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
		 *   NOTE: The fast order parameter is calculated as S2NH/S2s\n 
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
		 * EGAF parameters\n 
		 *   [0] tau slow\n 
		 *   [1] tau fast\n 
		 *   [2-4] alpha, beta, gamma deflections for slow motions\n 
		 *   [5] fast motion order parameter\n
		 * EGAFT parameters\n 
		 *   [0] tau slow\n 
		 *   [1] tau fast\n 
		 *   [2-4] alpha, beta, gamma deflections for slow motions\n 
		 *   [5] order parameter for fast motions\n 
		 *   [6] activation energy for slow motion\n 
		 *   [7] activation energy for fast motion\n 
		 */
		if (model == MOD_SMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
		} else if (model == MOD_EMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = resid->S2NH + (1 - resid->S2NH)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
			//printf("RUN: %f, %Le, %Le, %Le\n", resid->S2NH, opts[0], opts[1], opts[2]);
		} else if (model == MOD_EMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
			opts[1] = resid->S2NH + (1 - resid->S2NH)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
			opts[3] = (rand()%60000)/1.;
			opts[4] = (rand()%60000)/1.;
			//printf("RUN: %f, %Le, %Le, %Le\n", resid->S2NH, opts[0], opts[1], opts[2]);
		} else if (model == MOD_DEMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = resid->S2NH + (1 - resid->S2NH)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
			opts[3] = resid->S2NH + (1 - resid->S2NH)*((rand() % 100) / 100.);
		} else if (model == MOD_DEMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
			opts[1] = resid->S2NH + (1 - resid->S2NH)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
			opts[3] = resid->S2NH + (1 - resid->S2NH)*((rand() % 100) / 100.);
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
		} else if (model == MOD_EGAF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
			for (k = 2; k <= 4; k++) {
				// 15 degrees = 0.26180 radians
				opts[k] = ((rand () % 250)/1000.);
				//printf("%d %Lf\n", k, opts[k] * (180. / M_PI));
			}
			opts[5] = resid->S2NH + (1 - resid->S2NH)*((rand() % 100) / 100.);
			
		} else if (model == MOD_EGAFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
			opts[1] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
			for (k = 2; k <= 4; k++) {
				// 15 degrees = 0.26180 radians
				opts[k] = ((rand () % 250)/1000.);
				//printf("%d %Lf\n", k, opts[k] * (180. / M_PI));
			}
			opts[5] = resid->S2NH + (1 - resid->S2NH)*((rand() % 100) / 100.);
			opts[6] = (rand()%60000)/1.;
			opts[7] = (rand()%60000)/1.;
		}

		if (m->or_variation == VARIANT_A && m->rdc == RDC_ON) {
			/* eg for MOD_GAF, params = 8 + 3. opts[7] is full, so we want to put
			 * alpha in opts[params-3] and beta in opts[params-2] and gamma in opts[params-1];
			 */
			opts[params - 6] = 0; // papbC
			opts[params - 5] = 0; // papbN
			opts[params - 4] = (rand() % 20000) / 1.; // kex
			opts[params - 3] = 0; // alpha
			opts[params - 2] = 0; // beta
			opts[params - 1] = 0; // gamma
		} else if (m->or_variation == VARIANT_A) {
			opts[params - 3] = 0; // alpha
			opts[params - 2] = 0; // beta
			opts[params - 1] = 0; // gamma
		} else if (m->rdc == RDC_ON) {
			opts[params - 3] = 0; // papbC
			opts[params - 2] = 0; // papbN
			opts[params - 1] = (rand() % 20000) / 1.; // kex
		}
		
		//printf("%Le, %Le\n", opts[0], opts[1]);
		double val = simplex(optimize_chisq, opts, 1.0e-16, 1, resid, m);
		if (val >= 1000000 || val < 0) {
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
	start_time = time(0);
	srand((unsigned int)time(NULL));
	initialise_dwig(HALF_PI, Dwig);


	char system_file[255] = "";
	unsigned int i;
	int err_mod = 0;
	for (i = 1; i < (unsigned int) argc; i++) {
		//printf("%s\n", argv[i]);
		if (strcmp(argv[i], "-e") == 0)
			err_mod = 1;
		else if (strcmp(argv[i], "verify") == 0)
			verify_all();
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
	if (m.error_mode == 1 && m.n_error_iter == 0) {
		printf("Please provide number of error iterations\n");
		ret = -1;
	}

	if (m.model == MOD_UNDEFINED) {
		printf("Model undefined\n");
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

	unsigned int current_residue = 0;
	unsigned int n_spawns = 0;
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
	unsigned int l = 0;
	for (l = 0; l < n_spawns; l++) {
		for (i=0; i<m.nthreads; ++i) {
			RRA[i].i = current_residue + i;
			if (current_residue + i >= m.n_residues)
				continue;
			RRA[i].resid = &(m.residues[current_residue + i]);
			RRA[i].model = &m;
			RRA[i].n_iter = m.n_iter;
			strcpy(RRA[i].outputdir, m.outputdir);
			//printf("spawning thread %d (residue %d)\n", i, current_residue + i);
			rc = pthread_create(&threads[i], &threadattr, run_residue, (void *) &RRA[i]);
			if (rc != 0) {
				ERROR("Failed to spawn thread %d. Crashing...", current_residue+i);
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
		ERROR("%s not found.", filename);
		free_all(&m);
		return -1;
	}
	unsigned int params = m.params;
	FILE * ep = NULL;
	
	
	
	if (m.error_mode == 1) {
		sprintf(filename, "%s/errors.dat", m.outputdir);
		ep = fopen(filename, "w");
		if (ep == NULL) {
			ERROR("%s not found.", filename);
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
				RRA[i].model = &m;
				RRA[i].n_iter = m.n_error_iter;
				strcpy(RRA[i].outputdir, m.outputdir);
				//strcpy(RRA[i].outputdir, m.outputdir);
				//printf("spawning thread %d (residue %d)\n", i, current_residue + i);
				rc = pthread_create(&threads[i], &threadattr, calc_errors, (void *) &RRA[i]);
				if (rc != 0) {
					ERROR("Failed to spawn thread %d. Crashing...", current_residue+i);
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

	unsigned int k;
	FILE * gaf;
	if (m.model == MOD_GAF || m.model == MOD_GAFT || m.model == MOD_EGAF || m.model == MOD_EGAFT) {
		sprintf(filename, "%s/gaf.dat", m.outputdir);
		gaf = fopen(filename, "w");
		if (gaf == NULL) {
			ERROR("%s not found.", filename);
			free_all(&m);
			return -1;
		}
	}

	FILE * orderparams;
	sprintf(filename, "%s/orderparams.dat", m.outputdir);
	orderparams = fopen(filename, "w");
	if (orderparams == NULL) {
		ERROR("%s not found.", filename);
		free_all(&m);
		return -1;
	}

	printf("Outputting Files...\n");
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
		fprintf(fp, "%d\t%f\t%f", l+1, m.residues[l].S2NH, m.residues[l].min_val);
		if (m.error_mode == 1) {
			fprintf(ep, "%d\t%f\t%f", l+1, m.residues[l].S2NH, m.residues[l].min_val);
		}
		for (i = 0; i < params; i++) {
			fprintf(fp, "\t%Le", m.residues[l].parameters[i]);
			if (m.error_mode == 1) {
				//printf("%Le\n", m.residues[l].errors_std[i]);
				fprintf(ep, "\t%Le\t%Le", m.residues[l].parameters[i], 2 * m.residues[l].errors_std[i]);
			}
			/* WARNING: I'm printing the actual minimized parameters with the errors from calculation.
			 * The error_means are generally not minimal. */
		}

		fprintf(fp, "\t%d\n", m.residues[l].error_calcs);
		if (m.error_mode == 1) {
			fprintf(ep, "\n");
		}
		sprintf(file, "%s/backcalc_%d.dat", m.outputdir, l+1);
		back_calculate((m.residues[l].parameters), &(m.residues[l]), &m, file, params);
		
		double alpha, beta, gamma;
		if (m.or_variation == VARIANT_A) {
			alpha = (double) m.residues[l].parameters[m.params-3];
			beta = (double) m.residues[l].parameters[m.params-2];
			gamma = (double) m.residues[l].parameters[m.params-1];
			for (i = 0; i < N_OR; i++) {
				calculate_Y2(&(m.residues[l].orients[i]));
				rotate_Y2(&(m.residues[l].orients[i]), alpha, beta, gamma);
			}
		}

		if ((m.model == MOD_GAF || m.model == MOD_GAFT) && gaf != NULL) {
			double S2NHs, S2NHf, S2CHs, S2CHf, S2CNs, S2CNf; 
			long double sigs[3] = {m.residues[l].parameters[2], m.residues[l].parameters[3], m.residues[l].parameters[4]};
			long double sigf[3] = {m.residues[l].parameters[5], m.residues[l].parameters[6], m.residues[l].parameters[7]};
			struct Orient *Os[] = {&(m.residues[l].orients[OR_NH]), &(m.residues[l].orients[OR_CNH]), &(m.residues[l].orients[OR_CN])};
			double *S2s[] = {&S2NHs, &S2CHs, &S2CNs};
			double *S2f[] = {&S2NHf, &S2CHf, &S2CNf};
			
			printf("%f, %f :: %f\n", S2CHs, S2CHf, S2CHs * S2CHf);
			GAF_S2(sigs, Os, Os, S2s, 3, MODE_REAL);
			GAF_S2(sigf, Os, Os, S2f, 3, MODE_REAL);
			fprintf(orderparams, "%d\t", l+1);
			fprintf(orderparams, "%f\t%f\t%f\t", S2NHs*S2NHf, m.residues[l].S2NH, m.residues[l].S2NHe);
			fprintf(orderparams, "%f\t%f\t%f\t", S2CHs*S2CHf, m.residues[l].S2CH, m.residues[l].S2CHe);
			fprintf(orderparams, "%f\t%f\t%f\n", S2CNs*S2CNf, m.residues[l].S2CN, m.residues[l].S2CNe);

			long double taus = m.residues[l].parameters[0];
			long double tauf = m.residues[l].parameters[1];
			if (m.model == MOD_GAFT) {
				long double Eas = m.residues[l].parameters[8];
				long double Eaf = m.residues[l].parameters[9];
				taus *= expl(Eas / (RYD * 300));
				tauf *= expl(Eaf / (RYD * 300));
			}

			fprintf(gaf, "%d\t%Le\t%f\t%Le\t%f\n", l+1, taus, S2NHs, tauf, S2NHf);
		} else if ((m.model == MOD_EGAF || m.model == MOD_EGAFT) && gaf != NULL) {
			double S2NHs, S2NHf = (double) m.residues[l].parameters[5];
			double *S2[] = {&S2NHs};
			/* Approximate as just the S2NHs and S2NHf */
			struct Orient *As[] = {&(m.residues[l].orients[OR_NH])};
			long double sigs[3] = {m.residues[l].parameters[2], m.residues[l].parameters[3], m.residues[l].parameters[4]};
			GAF_S2(sigs, As, As, S2, 1, MODE_REAL);


			fprintf(orderparams, "%d\t", l+1);
			fprintf(orderparams, "%f\t%f\t%f\t", S2NHs*S2NHf, m.residues[l].S2NH, m.residues[l].S2NHe);
			fprintf(orderparams, "%f\t%f\t%f\t", S2NHs*S2NHf, m.residues[l].S2NH, m.residues[l].S2NHe);
			fprintf(orderparams, "%f\t%f\t%f\n", S2NHs*S2NHf, m.residues[l].S2NH, m.residues[l].S2NHe);


			long double taus = m.residues[l].parameters[0];
			long double tauf = m.residues[l].parameters[1];
			if (m.model == MOD_EGAFT) {
				long double Eas = m.residues[l].parameters[6];
				long double Eaf = m.residues[l].parameters[7];
				taus *= expl(Eas / (RYD * 300));
				tauf *= expl(Eaf / (RYD * 300));
			}
			fprintf(gaf, "%d\t%Le\t%f\t%Le\t%f\n", l+1, taus, S2NHs, tauf, S2NHf);
			/* This gives rise to an uninitialized warning, but this is initialized on line 445 inside another if so it should be fine. */
		} else if (m.model == MOD_DEMF || m.model == MOD_DEMFT) {
			double S2NHs, S2NHf;
			S2NHs = (double) m.residues[l].parameters[1];
			S2NHf = (double) m.residues[l].parameters[3];
			fprintf(orderparams, "%d\t", l+1);
			fprintf(orderparams, "%f\t%f\t%f\t", S2NHs*S2NHf, m.residues[l].S2NH, m.residues[l].S2NHe);
			fprintf(orderparams, "%f\t%f\t%f\t", S2NHs*S2NHf, m.residues[l].S2NH, m.residues[l].S2NHe);
			fprintf(orderparams, "%f\t%f\t%f\n", S2NHs*S2NHf, m.residues[l].S2NH, m.residues[l].S2NHe);
		} else if (m.model == MOD_SMF || m.model == MOD_SMFT) {
			double S2 = (double) m.residues[l].parameters[1];

			fprintf(orderparams, "%d\t", l+1);
			fprintf(orderparams, "%f\t%f\t%f\t", S2, m.residues[l].S2NH, m.residues[l].S2NHe);
			fprintf(orderparams, "%f\t%f\t%f\t", S2, m.residues[l].S2NH, m.residues[l].S2NHe);
			fprintf(orderparams, "%f\t%f\t%f\n", S2, m.residues[l].S2NH, m.residues[l].S2NHe);
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
