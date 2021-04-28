/*
 * Copyright (c) 2021 Ben Tatman, University of Warwick
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the Sosftware
 * is furnished to do so, subject to the following conditions:
 *
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software
 * - Any academic work deriving from the Software shall cite [CITATION].
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUR OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * runners.c
 *
 * Provides functions to run models and output errors on a per residue basis.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "datatypes.h"
#include "anneal.h"
#include "chisq.h"
#include "crosen.h"
#include "errors.h"
#include "runners.h"
#include "models/model.h"

/**
 * Operates residue optimization. Generates random parameter guesses and passes these to the simplex function.
 * @param input
 *  Pointer to rrarg containing thread information
 */
int run_residue(struct Model *m, unsigned int residue) {
	struct Residue * resid = &(m->residues[residue]);



	char outputdir[255];
	strcpy(outputdir, m->outputdir);
	FILE *fp;
	char filename[300];

	sprintf(filename, "%s/residue_%d.dat", outputdir, residue+1);
	fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}


	//printf("RESIDUE %d\n", i+1);
	//printf("Number of relaxations: %d\n", resid->n_relaxation);
	unsigned int l, k, q;
	unsigned int params = m->params;
	//printf("%d\n", params);
	Decimal *opts;
	opts = (Decimal *) malloc (sizeof(Decimal) * params);
	Decimal *anneal_pars = (Decimal *) malloc(sizeof(Decimal) * params);
	resid->parameters = (Decimal *) malloc (sizeof(Decimal) * params);
	if (resid->parameters == NULL || opts == NULL) {
		printf("Allocation failed.\n");
	}
	resid->min_val = MIN_VAL;

	Decimal *minv = (Decimal *) malloc(sizeof(Decimal) * params);
	Decimal *maxv = (Decimal *) malloc(sizeof(Decimal) * params);

	setup_paramlims(m, resid->S2NH, minv, maxv);
	Decimal val;

	if (resid->ignore == 1) {
		free(minv);
		free(maxv);
		free(anneal_pars);
		free(opts);
		fprintf(fp, "%d, %lf, -1, -1\n", residue + 1, 1000.);
		fclose(fp);
		return -1;
	}

	for (l = 0; l < m->n_anneal_iter; l++) {
		val = anneal(optimize_chisq, anneal_pars, minv, maxv, params, m->anneal_temp, 1000, m->anneal_wobb, m->anneal_therm, m->anneal_restart, NULL, MODE_RANDOM_START +MODE_RANDOM_RESTART, resid, m);
		for (q = 0; q < m->n_nm_iter; q++) {
			for (k = 0; k < params; k++) {
				opts[k] = anneal_pars[k];
			}
			val = simplex(optimize_chisq, opts, 1.0e-16, 1, resid, m);
			if (val >= 1000000 || val < 0) {
				val = -1;
				for (k = 0; k < params; k++) {
					opts[k] = -1;
				}
			}
			fprintf(fp, "%d\t%lf", residue + 1, val);
			for (k = 0; k < params; k++) {
				fprintf(fp, "\t%le", opts[k]);
			}

			if (val < resid->min_val && val != -1) {
				//printf("New lowest %lf\n", val);
				resid->min_val = val;
				for (k = 0; k < params; k++) {
					resid->parameters[k] = opts[k];
				}
			}
			fprintf(fp, "\n");
		}
	}
	free(minv);
	free(maxv);
	free(anneal_pars);
	free(opts);
	fclose(fp);
	return 0;
}

int run_fitting(struct Model *m) {
	unsigned int i;
	for (i = m->proc_start; i < m->proc_end; i++) {
		printf("\tWorker %d/%d: Residue %d.\n", m->myid + 1, m->numprocs, i + 1);
		run_residue(m, i);
	}
	return 1;
}

int run_errors(struct Model *m) {
	unsigned int i;
	for (i = m->proc_start; i < m->proc_end; i++) {
		printf("\tWorker %d/%d: Residue %d (errors).\n", m->myid + 1, m->numprocs, i+1);
		calc_errors(m, i);
	}
	return 1;
}

int print_errors(struct Model *m) {
	FILE * ep = NULL;
	char filename[300];
	unsigned int l, k, i;
	sprintf(filename, "%s/errors_proc%d.dat", m->outputdir, m->myid);
	ep = fopen(filename, "w");
	if (ep == NULL) {
		ERROR("%s not found.", filename);
		free_all(m);
		return -1;
	}
	for (l = m->proc_start; l < m->proc_end; l++) {
		for (k = 0; k < m->params; k++) {
			calc_statistics(m->residues[l].error_params[k], m->residues[l].error_calcs, &(m->residues[l].errors_mean[k]), &(m->residues[l].errors_std[k]));
		}

		fprintf(ep, "%d\t%lf\t%lf", l+1, m->residues[l].S2NH, m->residues[l].min_val);

		for (i = 0; i < m->params; i++) {

			fprintf(ep, "\t%le\t%le", m->residues[l].parameters[i], 2 * m->residues[l].errors_std[i]);
			/* WARNING: I'm printing the actual minimized parameters with the errors from calculation.
			 * The error_means are generally not minimal. */
		}

		fprintf(ep, "\n");
	}
	fclose(ep);
	return 1;
}

int print_residues(struct Model *m) {
	FILE * fp;
	char filename[300];
	unsigned int l, i;
	sprintf(filename, "%s/final_proc%d.dat", m->outputdir, m->myid);
	fp = fopen(filename, "w");
	if (fp == NULL) {
		ERROR("%s not found.", filename);
		free_all(m);
		return -1;
	}
	for (l = m->proc_start; l < m->proc_end; l++) {
		if (m->residues[l].min_val == MIN_VAL) {
			m->residues[l].min_val = -1.;
			for (i = 0; i < m->params; i++) {
				m->residues[l].parameters[i] = -1.;
			}
		}
		fprintf(fp, "%d\t%lf\t%lf", l+1, m->residues[l].S2NH, m->residues[l].min_val);
		for (i = 0; i < m->params; i++) {
			fprintf(fp, "\t%le", m->residues[l].parameters[i]);
		}
		fprintf(fp, "\t%d\n", m->residues[l].error_calcs);
	}
	fclose(fp);
	return 1;

}

int print_gaf(struct Model *m) {
	FILE * gaf;
	unsigned int l, i;
	char filename[300];
	sprintf(filename, "%s/gaf_proc%d.dat", m->outputdir, m->myid);
	gaf = fopen(filename, "w");
	if (gaf == NULL) {
		ERROR("%s not found.", filename);
		free_all(m);
		return -1;
	}
	FILE * orderparams;
	sprintf(filename, "%s/orderparams_proc%d.dat", m->outputdir, m->myid);
	orderparams = fopen(filename, "w");
	if (orderparams == NULL) {
		ERROR("%s not found.", filename);
		free_all(m);
		return -1;
	}
	int ignore = -1;
	for (l = m->proc_start; l < m->proc_end; l++) {
		Decimal alpha, beta, gamma;
		if (m->or_variation == VARIANT_A) {
			alpha = (Decimal) m->residues[l].parameters[m->params - 3];
			beta = (Decimal) m->residues[l].parameters[m->params - 2];
			gamma = (Decimal) m->residues[l].parameters[m->params - 1];
			for (i = 0; i < N_OR; i++) {
				calculate_Y2(&(m->residues[l].orients[i]));
				rotate_Y2(&(m->residues[l].orients[i]), alpha, beta, gamma);
			}
		}

		struct BCParameters pars;
		int otb = opts_to_bcpars(m->residues[l].parameters, &pars, m, &(m->residues[l]), &ignore);
		if (otb != 0)
			return -1;

		double S2NH, S2CH, S2CC, S2CN;
		S2NH = pars.S2NHs * pars.S2NHf * pars.S2uf;
		S2CH = pars.S2CHs * pars.S2CHf * pars.S2uf;
		S2CC = pars.S2CCAps * pars.S2CCApf * pars.S2uf;
		S2CN = pars.S2CNs * pars.S2CNf * pars.S2uf;

		fprintf(orderparams, "%d\t", l + 1);
		fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NH, m->residues[l].S2NH, m->residues[l].S2NHe);
		fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CH, m->residues[l].S2CH, m->residues[l].S2CHe);
		fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CN, m->residues[l].S2CN, m->residues[l].S2CNe);
		fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CC, m->residues[l].S2CC, m->residues[l].S2CCe);
		fprintf(orderparams, "\n");

		fprintf(gaf, "%d\t%le\t%lf\t%le\t%lf\n", l + 1, pars.taus, pars.S2NHs, pars.tauf, pars.S2NHf);
	}
	fclose(gaf);
	fclose(orderparams);
	return 1;
}

int print_backcalcs(struct Model *m) {
	char filename[300];
	unsigned int l;
	for (l = m->proc_start; l < m->proc_end; l++) {
		sprintf(filename, "%s/backcalc_%d.dat", m->outputdir, l+1);
		back_calculate((m->residues[l].parameters), &(m->residues[l]), m, filename, m->params);
	}
	return 1;
}

