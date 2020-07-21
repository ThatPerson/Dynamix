#include <stdio.h>
#include "datatypes.c"
#include "read_data.c"
#include "models.c"
#include "chisq.c"
#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
#include <time.h>
#include <pthread.h>



void * run_residue(void *input) {
	int i = ((struct rrargs*)input)->i;
	struct Residue * resid = ((struct rrargs*)input)->resid;
	int model = ((struct rrargs*)input)->model;
	int n_iter = ((struct rrargs*)input)->n_iter;
	char outputdir[255];
	strcpy(outputdir, ((struct rrargs*)input)->outputdir);
	FILE *fp;
	char line[255];
	size_t len=255;
	char filename[255];
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
}




int main(void) {
	srand(time(NULL));
	struct Model m;
	int ret = read_system_file("system/model.dx", &m);
	printf("%d\n", ret);
	//print_system(&m);
	
	int i;
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
	char filename[255];
	FILE * fp;
	sprintf(filename, "%s/final.dat", m.outputdir);
	fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	int params;
	switch (m.model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT:params = 5; break;
		case MOD_SMFT:params = 3; break;
		default: params = 0; break;
	}
	long double c = -1;
	char file[255];
	for (l = 0; l < m.n_residues; l++) {
		if (m.residues[l].min_val == MIN_VAL) {
			fprintf(fp, "%d\t%f\t%f", l + 1, m.residues[l].S2_dipolar, -1.);
			for (i = 0; i < params; i++) {
				fprintf(fp, "\t%Le", c);
			}
		} else {
			fprintf(fp, "%d\t%f\t%f", l, m.residues[l].S2_dipolar, m.residues[l].min_val);
			for (i = 0; i < params; i++) {
				fprintf(fp, "\t%Le", m.residues[l].parameters[i]);
			}
		}
		fprintf(fp, "\n");
		sprintf(file, "%s/backcalc_%d.dat", m.outputdir, l);
		back_calculate((m.residues[l].parameters), &(m.residues[l]), m.model, file);
	}
	free_all(&m);
	return 1;
}
