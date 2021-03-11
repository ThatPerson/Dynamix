
Decimal optimize_global_chisq(Decimal * opts, struct Residue *r, struct Model * m, unsigned int params) {
	/* opts is a pointer to an array containing;
	 *
	 *  for SMF, [tau, S2]
	 */
    (void) r;
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
	unsigned int params = m->params;
	//printf("%d\n", params);
	Decimal *opts;
	opts = (Decimal *) malloc (sizeof(Decimal) * params);
    Decimal *anneal_pars = (Decimal *) malloc(sizeof(Decimal) * params);
    Decimal *minv = (Decimal *) malloc(sizeof(Decimal) * params);
    Decimal *maxv = (Decimal *) malloc(sizeof(Decimal) * params);

    setup_paramlims(m, resid->S2NH, minv, maxv);

    if (resid->ignore == 1) {
        free(minv);
        free(maxv);
        free(anneal_pars);
        free(opts);
        fprintf(fp, "%d, %lf, -1, -1\n", i + 1, 1000.);
        fclose(fp);
        return -1;
    }

    unsigned int q;
    Decimal val;
    anneal(optimize_global_chisq, anneal_pars, minv, maxv, m->params, m->anneal_temp, 1000, m->anneal_wobb, m->anneal_therm, m->anneal_restart, NULL, MODE_RANDOM_RESTART+MODE_RANDOM_START, resid, m);

	//printf("%le, %le\n", opts[0], opts[1]);
	for (q = 0; q < m->n_nm_iter; q++) {
	    for (k = 0; k < m->params; k++) {
	        opts[k] = anneal_pars[k];
	    }
        val = simplex(optimize_global_chisq, opts, 1.0e-16, 1, resid, m);
        if (val >= 1000000 || val < 0) {
            val = -1;
            for (k = 0; k < m->params; k++) {
                opts[k] = -1;
            }
        }

        fprintf(fp, "%d\t%lf", i + 1, val);
        for (k = 0; k < params; k++) {
            fprintf(fp, "\t%le", opts[k]);
        }
        fprintf(fp, "\n");
        unsigned int j;
        for (j = 0; j < m->n_residues; j++) {
            if (val < m->residues[j].min_val && val != -1) {
                //printf("New lowest %lf\n", val);
                m->residues[j].min_val = val;
                for (k = 0; k < params; k++) {
                    m->residues[j].parameters[k] = opts[k];
                }
            }
        }
    }
	free(minv);
	free(maxv);


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
	Decimal temp_R;
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
	for (l = 0; l < m->n_anneal_iter; l++) {
		printf("\tRunning iteration %d\n", l);
		run_global_iteration(m, l);
	}
	return 1;
}

