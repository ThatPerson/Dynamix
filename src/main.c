#include <stdio.h>
#include "datatypes.c"
#include "read_data.c"
#include "models.c"
#include "chisq.c"
#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
#include "errors.c"
#include <time.h>
#include <pthread.h>



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
		if (model == MOD_SMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8);
			opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
		} else if (model == MOD_EMF) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -8);			
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from 0.5 to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -11);
			//printf("RUN: %f, %Le, %Le, %Le\n", resid->S2_dipolar, opts[0], opts[1], opts[2]);
		} else if (model == MOD_EMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15);			
			opts[1] = resid->S2_dipolar + (1 - resid->S2_dipolar)*((rand() % 100) / 100.); // random number from 0.5 to 1
			opts[2] = ((rand() % 100)/100.) * powl(10, -20);
			opts[3] = (rand()%60000)/1.;
			opts[4] = (rand()%60000)/1.;
			//printf("RUN: %f, %Le, %Le, %Le\n", resid->S2_dipolar, opts[0], opts[1], opts[2]);
		} else if (model == MOD_SMFT) {
			opts[0] = ((rand() % 100)/100.) * powl(10, -15);
			opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
			opts[2] = (rand()%60000)/1.;
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




int main(int argc, char * argv[]) {
	srand(time(NULL));
	int error_mode = 0;
	char system_file[255] = "";
	int i;
	for (i = 1; i < argc; i++) {
		//printf("%s\n", argv[i]);
		if (strcmp(argv[i], "-e") == 0)
			error_mode = 1;
		else
			strcpy(system_file, argv[i]);
	}
	
	if (strcmp(system_file, "") == 0){ 
		printf("Please provide system file.\n");
		exit(-1);
	}
	
	
	
	struct Model m;
	int ret = read_system_file(system_file, &m);
	printf("%d\n", ret);
	//print_system(&m);
	
	/*for (i = 0; i < 10; i++) {
		long double taus = (long double) ((rand() % 100) / 100.) * pow(10, -(5+(rand()%5)));
		long double S2s = (long double) ((rand()%100)/100.);
		long double tauf = (long double) ((rand() % 100) / 100.) * pow(10, -(7 + (rand()%5)));
		printf("%f\n", EMF_15NR2(&(m.residues[40]), &(m.residues[40].relaxation[19]), taus, S2s, tauf));
		printf("\n\n\n\n");
	}
	
	exit(-1);*/
	
	//int i;
	
	
	pthread_t threads[NTHREADS];
	int rc;
	
	int current_residue = 0;
	int n_spawns = 0;
	/* if we have 56 residues and 4 threads then we need 
	 * 56 / 4 spawn events (= 14). Add 1 in case (eg for 57).
	 * Then loop over threads, increment current_residue and assign pointers.
	 */
	n_spawns = m.n_residues / NTHREADS;
	n_spawns ++;
	//printf("%d spawns\n", n_spawns);

	struct rrargs RRA[NTHREADS]; // = (struct rrargs *)malloc(sizeof(struct rrargs));
	printf("Optimizing...\n");
	/* spawn the threads */
	int l = 0;
	for (l = 0; l < n_spawns; l++) {
		for (i=0; i<NTHREADS; ++i) {
			RRA[i].i = current_residue + i;
			if (current_residue + i >= m.n_residues)
				continue;
			RRA[i].resid = &(m.residues[current_residue + i]);
			RRA[i].model = m.model;
			RRA[i].n_iter = m.n_iter;
			strcpy(RRA[i].outputdir, m.outputdir);
			//printf("spawning thread %d (residue %d)\n", i, current_residue + i);
			rc = pthread_create(&threads[i], NULL, run_residue, (void *) &RRA[i]);
		}
		current_residue += NTHREADS;

	
		for (i=0; i<NTHREADS; ++i) {
			rc = pthread_join(threads[i], NULL);
		}
	}
	char filename[300];
	FILE * fp;
	sprintf(filename, "%s/final.dat", m.outputdir);
	fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	int params;
	FILE * ep = NULL;
	
	switch (m.model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT:params = 5; break;
		case MOD_SMFT:params = 3; break;
		default: params = 0; break;
	}
	
	if (error_mode == 1) {
		sprintf(filename, "%s/errors.dat", m.outputdir);
		ep = fopen(filename, "w");
		if (ep == NULL) {
			printf("%s not found.\n", filename);
			return -1;
		}
	
		
		
		printf("Calculating Errors...\n");
		current_residue = 0;
		for (l = 0; l < n_spawns; l++) {
			for (i=0; i<NTHREADS; ++i) {
				RRA[i].i = current_residue + i;
				if (current_residue + i >= m.n_residues)
					continue;
				RRA[i].resid = &(m.residues[current_residue + i]);
				RRA[i].model = m.model;
				RRA[i].n_iter = m.n_error_iter;
				//strcpy(RRA[i].outputdir, m.outputdir);
				//printf("spawning thread %d (residue %d)\n", i, current_residue + i);
				rc = pthread_create(&threads[i], NULL, calc_errors, (void *) &RRA[i]);
			}
			current_residue += NTHREADS;

		
			for (i=0; i<NTHREADS; ++i) {
				rc = pthread_join(threads[i], NULL);
			}
		}
	}
	
	int k;
	
	
	printf("Outputting Files...\n");
	long double c = -1;
	char file[300];
	for (l = 0; l < m.n_residues; l++) {
		if (error_mode == 1) {
			for (k = 0; k < params; k++) {
				calc_statistics(m.residues[l].error_params[k], m.n_error_iter, &(m.residues[l].errors_mean[k]), &(m.residues[l].errors_std[k])); 
			}
		}
		if (m.residues[l].min_val == MIN_VAL) {
			m.residues[l].min_val = -1.;
			for (i = 0; i < params; i++) {
				m.residues[l].parameters[i] = -1.;
			}
		}
		fprintf(fp, "%d\t%f\t%f", l, m.residues[l].S2_dipolar, m.residues[l].min_val);
		if (error_mode == 1) {
			fprintf(ep, "%d, %f, %f", l, m.residues[l].S2_dipolar, m.residues[l].min_val);
		}
		for (i = 0; i < params; i++) {
			fprintf(fp, "\t%Le", m.residues[l].parameters[i]);
			if (error_mode == 1) {
				fprintf(ep, ", %Le, %Le", m.residues[l].parameters[i], 2 * m.residues[l].errors_std[i]);
			}
			/* WARNING: I'm printing the actual minimized parameters with the errors from calculation.
			 * The error_means are generally not minimal. */
		}
		
		fprintf(fp, "\n");
		if (error_mode == 1) {
			fprintf(ep, "\n");
		}
		sprintf(file, "%s/backcalc_%d.dat", m.outputdir, l+1);
		back_calculate((m.residues[l].parameters), &(m.residues[l]), m.model, file);
	}
	
	
	fclose(fp);
	if (error_mode == 1)
		fclose(ep);
	free_all(&m);
	return 1;
}
