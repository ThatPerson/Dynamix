
double optimize_global_chisq(long double * opts, struct Residue *r, struct Model * m, unsigned int params) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */

	unsigned int l, c = 0;
	double chisq = 0;
	for (l = 0; l < m->n_residues; l++) {
	//	double optimize_chisq(long double * opts, struct Residue * resid, struct Model * m, unsigned int params) {
		if (m->residues[l].ignore == 1)
			continue;
		chisq += optimize_chisq(opts, &(m->residues[l]), m, params);
		c++;
	}
	
	return (double) (chisq / c); 
	/* normalise to number of relaxation measurements - otherwise when using like 85 the chisq becomes huge which hinders convergence */
}

pthread_mutex_t write_lock;
pthread_mutex_t param_lock;

void *run_global_iteration(void *input) {
	unsigned int i = ((struct rrargs*)input)->i;
	printf("\tThread %d alive...\n", i+1);
	struct Model *m = ((struct rrargs*)input)->model;
	struct Residue *resid = NULL;
	FILE *fp;
	char filename[300];
	sprintf(filename, "%s/fitting.dat", m->outputdir);
	fp = fopen(filename, "a");
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
	
	

		/*if (resid->ignore == 1) {
			fprintf(fp, "%d, %f, -1, -1\n", i+1, 1000.);
			fclose(fp);
			free(opts);
			return NULL;
		}*/
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
	unsigned int model = m->model;
	if (model == MOD_SMF) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
		opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
	} else if (model == MOD_EMF) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
		opts[1] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
		opts[2] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
		//printf("RUN: %f, %Le, %Le, %Le\n", 0.7, opts[0], opts[1], opts[2]);
	} else if (model == MOD_EMFT) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
		opts[1] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
		opts[2] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
		opts[3] = (rand()%60000)/1.;
		opts[4] = (rand()%60000)/1.;
		//printf("RUN: %f, %Le, %Le, %Le\n", 0.7, opts[0], opts[1], opts[2]);
	} else if (model == MOD_DEMF) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
		opts[1] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
		opts[2] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
		opts[3] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.);
	} else if (model == MOD_DEMFT) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
		opts[1] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
		opts[2] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
		opts[3] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.);
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
	} else if (model == MOD_AIMF) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
		opts[1] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
		for (k = 2; k <= 7; k++) {
			opts[k] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.);
			//printf("%d %Lf\n", k, opts[k] * (180. / M_PI));
		}
	} else if (model == MOD_AIMFT) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
		opts[1] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
		for (k = 2; k <= 7; k++) {
			opts[k] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.);
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
		opts[5] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.);
		
	} else if (model == MOD_EGAFT) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
		opts[1] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
		for (k = 2; k <= 4; k++) {
			// 15 degrees = 0.26180 radians
			opts[k] = ((rand () % 250)/1000.);
			//printf("%d %Lf\n", k, opts[k] * (180. / M_PI));
		}
		opts[5] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.);
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
	double val = simplex(optimize_global_chisq, opts, 1.0e-16, 1, resid, m);
	if (val >= 1000000 || val < 0) {
		val = -1;
		for (k = 0; k < m->params; k++) {
			opts[k] = -1;
		}
	}
	pthread_mutex_lock(&write_lock);
	fprintf(fp, "%d\t%f", l+1, val);
	for (k = 0; k < params; k++) {
		fprintf(fp, "\t%Le", opts[k]);
	}
	fprintf(fp, "\n");
	pthread_mutex_unlock(&write_lock);
	
	pthread_mutex_lock(&param_lock);
	for (i = 0; i < m->n_residues; i++) {
		if (val < m->residues[i].min_val && val != -1) {
		//printf("New lowest %f\n", val);	
			m->residues[i].min_val = val;
			for (k = 0; k < params; k++) {
				m->residues[i].parameters[k] = opts[k];
			}
		}
	}
	pthread_mutex_unlock(&param_lock);


	free(opts);
	fclose(fp);
	return NULL;
}




/**
 * Operates residue optimization. Generates random parameter guesses and passes these to the simplex function.
 * @param input
 *  Pointer to rrarg containing thread information
 */
void * run_global(struct Model *m) {
	if ((pthread_mutex_init(&write_lock, NULL) != 0) || (pthread_mutex_init(&param_lock, NULL) != 0)) {
		ERROR("Mutex init failed.\n");
		free_all(m);
		exit(-1);
	}
	unsigned int l;
	for (l = 0; l < m->n_residues; l++) {
		m->residues[l].parameters = (long double *) malloc (sizeof(long double) * m->params);
		m->residues[l].min_val = MIN_VAL;
	}
	
	// launch
	pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * m->nthreads);
	pthread_attr_t threadattr;
	pthread_attr_init(&threadattr);
	pthread_attr_setstacksize(&threadattr, THREAD_STACK);
	int rc;

	unsigned int current_it = 0;
	unsigned int n_spawns = 0, i, current_residue;
	n_spawns = m->n_iter / m->nthreads;
	n_spawns ++;
	struct rrargs * RRA = (struct rrargs *)malloc(sizeof(struct rrargs) * m->nthreads);
	
	unsigned int params = m->params;


	current_residue = 0;
	for (l = 0; l < n_spawns; l++) {
		for (i=0; i<m->nthreads; ++i) {
			RRA[i].i = current_residue + i;
			if (current_residue + i >= m->n_residues)
				continue;
			RRA[i].resid = NULL;
			RRA[i].model = m;
			strcpy(RRA[i].outputdir, m->outputdir);
			//strcpy(RRA[i].outputdir, m.outputdir);
			//printf("spawning thread %d (residue %d)\n", i, current_residue + i);
			rc = pthread_create(&threads[i], &threadattr, run_global_iteration, (void *) &RRA[i]);
			if (rc != 0) {
				ERROR("Failed to spawn thread %d. Crashing...", current_residue+i);
				free_all(m);
				exit(-1);
			}
		}
		current_residue += m->nthreads;

		for (i=0; i<m->nthreads; ++i) {
			rc = pthread_join(threads[i], NULL);
		}
	}

	free(RRA);
	free(threads);
	
	pthread_mutex_destroy(&write_lock);
	pthread_mutex_destroy(&param_lock);
}

