
Decimal optimize_global_chisq(Decimal * opts, struct Residue *r, struct Model * m, unsigned int params) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */

	unsigned int l, c = 0;
	Decimal chisq = 0;
	for (l = 0; l < m->n_residues; l++) {
	//	Decimal optimize_chisq(Decimal * opts, struct Residue * resid, struct Model * m, unsigned int params) {
		if (m->residues[l].ignore == 1)
			continue;
		chisq += optimize_chisq(opts, &(m->residues[l]), m, params);
		c++;
	}
	
	return (Decimal) (chisq / c); 
	/* normalise to number of relaxation measurements - otherwise when using like 85 the chisq becomes huge which hinders convergence */
}

int run_global_iteration(struct Model *m, int i) {
	struct Residue *resid = NULL;
	FILE *fp;
	char filename[300];
	sprintf(filename, "%s/fitting.dat", m->outputdir);
	fp = fopen(filename, "a");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}

	//printf("RESIDUE %d\n", i+1);
	//printf("Number of relaxations: %d\n", resid->n_relaxation);
	unsigned int k;
	unsigned int params = 0;
	params = m->params;
	//printf("%d\n", params);
	Decimal *opts;
	opts = (Decimal *) malloc (sizeof(Decimal) * params);
	
	

		/*if (resid->ignore == 1) {
			fprintf(fp, "%d, %lf, -1, -1\n", i+1, 1000.);
			fclose(fp);
			free(opts);
			return NULL;
		}*/
		//Decimal opts[20] = {0, 0};
		
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
		//printf("RUN: %lf, %le, %le, %le\n", 0.7, opts[0], opts[1], opts[2]);
	} else if (model == MOD_EMFT) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -15 + T_S);
		opts[1] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.); // random number from s2 dipolar to 1
		opts[2] = ((rand() % 100)/100.) * powl(10, -20 + T_S);
		opts[3] = (rand()%60000)/1.;
		opts[4] = (rand()%60000)/1.;
		//printf("RUN: %lf, %le, %le, %le\n", 0.7, opts[0], opts[1], opts[2]);
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
			//printf("%d %lf\n", k, opts[k] * (180. / M_PI));
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
			//printf("%d %lf\n", k, opts[k] * (180. / M_PI));
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
			//printf("%d %lf\n", k, opts[k] * (180. / M_PI));
		}
		opts[5] = 0.7 + (1 - 0.7)*((rand() % 100) / 100.);
		
	} else if (model == MOD_EGAFT) {
		opts[0] = ((rand() % 100)/100.) * powl(10, -8 + T_S);
		opts[1] = ((rand() % 100)/100.) * powl(10, -11 + T_S);
		for (k = 2; k <= 4; k++) {
			// 15 degrees = 0.26180 radians
			opts[k] = ((rand () % 250)/1000.);
			//printf("%d %lf\n", k, opts[k] * (180. / M_PI));
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
	

	//printf("%le, %le\n", opts[0], opts[1]);
	Decimal val = simplex(optimize_global_chisq, opts, 1.0e-16, 1, resid, m);
	if (val >= 1000000 || val < 0) {
		val = -1;
		for (k = 0; k < m->params; k++) {
			opts[k] = -1;
		}
	}

	fprintf(fp, "%d\t%lf", i+1, val);
	for (k = 0; k < params; k++) {
		fprintf(fp, "\t%le", opts[k]);
	}
	fprintf(fp, "\n");

	for (i = 0; i < m->n_residues; i++) {
		if (val < m->residues[i].min_val && val != -1) {
		//printf("New lowest %lf\n", val);	
			m->residues[i].min_val = val;
			for (k = 0; k < params; k++) {
				m->residues[i].parameters[k] = opts[k];
			}
		}
	}



	free(opts);
	fclose(fp);
	return 1;
}

int calc_global_errors(struct Model *m) {
	FILE *errp;
	char filename[300];
	sprintf(filename, "%s/errors.dat", m->outputdir);
	errp = fopen(filename, "a");
	if (errp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	
	unsigned int l, k, i;
	struct Residue *resid;
	Decimal *opts;
	opts = (Decimal *) malloc(sizeof(Decimal) * m->params);
	
	
	int p = 0;
	Decimal temp_R = 0;
	int ignore = -1;
	
	for (i = 0; i < m->n_residues; i++) {
		resid = &(m->residues[i]);
		resid->S2NHb = resid->S2NH;
		resid->S2CHb = resid->S2CH;
		resid->S2CNb = resid->S2CN;
		resid->S2CCb = resid->S2CC;
		resid->errors_mean = (Decimal *) malloc (sizeof(Decimal) * m->params);
		resid->errors_std = (Decimal *) malloc(sizeof(Decimal) * m->params);
		resid->error_params = (Decimal **) malloc (sizeof(Decimal *) * m->params);
		for (k = 0; k < m->params; k++) {
			resid->error_params[k] = (Decimal *) malloc (sizeof(Decimal) * m->n_error_iter);
		}
	}
	

	
	for (l = 0; l < m->n_error_iter; l++) {
		printf("Iteration %d.\n", l);
		for (i = 0; i < m->n_residues; i++) {
			resid = &(m->residues[i]);
			resid->temp_relaxation = (resid->relaxation);
			resid->relaxation = NULL;
			resid->relaxation = (struct Relaxation *) malloc(sizeof(struct Relaxation) * resid->lim_relaxation);

			for (k = 0; k < resid->n_relaxation; k++) {
				resid->relaxation[k].field = resid->temp_relaxation[k].field;
				resid->relaxation[k].wr = resid->temp_relaxation[k].wr;
				resid->relaxation[k].w1 = resid->temp_relaxation[k].w1;
				resid->relaxation[k].type = resid->temp_relaxation[k].type;
				resid->relaxation[k].T = resid->temp_relaxation[k].T;
				
				temp_R = back_calc(resid->parameters, resid, &(resid->temp_relaxation[k]), m, &ignore);
				resid->relaxation[k].R = norm_rand(temp_R, (resid->temp_relaxation[k].Rerror/2.));
				LOG("%d %d %d prior: %lf, backcalc: %lf, new: %lf, error: %lf", i, l, k, resid->temp_relaxation[k].R, temp_R, resid->relaxation[k].R, resid->temp_relaxation[k].Rerror/2.);
				resid->relaxation[k].Rerror = resid->temp_relaxation[k].Rerror;
			}
			resid->S2NH = norm_rand(resid->S2NHb, (resid->S2NHe / 2.)); // divided by two to get 1 standard deviation
			resid->S2CH = norm_rand(resid->S2CHb, (resid->S2CHe / 2.));
			resid->S2CN = norm_rand(resid->S2CNb, (resid->S2CNe / 2.));
			resid->S2CC = norm_rand(resid->S2CCb, (resid->S2CCe / 2.));
		}
		for (k = 0; k < m->params; k++) {
			opts[k] = m->residues[0].parameters[k]; // doesn't matter which we take it from, all are the same.
		}

		resid = NULL;
		Decimal min = simplex(optimize_global_chisq, opts, 1.0e-16, 1, resid, m);
		LOG("%d %d min = %lf", i, l, min);
		fprintf(errp, "%d\t%lf", l+1, min);
		for (k = 0; k < m->params; k++) {
			for (i = 0; i < m->n_residues; i++)
				m->residues[i].error_params[k][p] = opts[k];
			fprintf(errp, "\t%le", opts[k]);
		}
		fprintf(errp, "\n");
		p++;
		for (i = 0; i < m->n_residues; i++) {
			resid = &(m->residues[i]);
			free(resid->relaxation);
			resid->relaxation = NULL;
			resid->relaxation = resid->temp_relaxation;
			resid->temp_relaxation = NULL;
		}
	}
	for (i = 0; i < m->n_residues; i++) {
		m->residues[i].S2NH = m->residues[i].S2NHb;
		m->residues[i].S2CH = m->residues[i].S2CHb;
		m->residues[i].S2CN = m->residues[i].S2CNb;
		m->residues[i].S2CC = m->residues[i].S2CCb;
		m->residues[i].error_calcs = p;
	}
	fclose(errp);
	free(opts);
	return 1;
}


/**
 * Operates residue optimization. Generates random parameter guesses and passes these to the simplex function.
 * @param input
 *  Pointer to rrarg containing thread information
 */
int run_global(struct Model *m) {
	if (m->numprocs > 1) {
		ERROR("Global GAF fitting _does not_ support MPI parallelisation at the moment.\n");
		free_all(m);
		exit(-1);
	}

	unsigned int l;
	for (l = 0; l < m->n_residues; l++) {
		m->residues[l].parameters = (Decimal *) malloc (sizeof(Decimal) * m->params);
		m->residues[l].min_val = MIN_VAL;
	}
	
	// launch
	for (l = 0; l < m->n_iter; l++) {
		printf("\tRunning iteration %d\n", l);
		run_global_iteration(m, l);
	}
	return 1;
}

