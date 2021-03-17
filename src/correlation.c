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
 * correlation.c
 *
 * Separate program from the main Dynamix. Reads in output files from Dynamix
 * and can then either perform backcalculations, calculate correlation functions
 * or spectral densities.
 */

#include <stdio.h>
#include <complex.h>
#include "datatypes.h"
#include "read_data.h"
#include <omp.h>
#include "models/model.h"
#include "chisq.h"
#include "crosen.h" // implementation of Nelder-Mead simplex algorithm
#include "errors.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include "verification.h"


/**
 * Reads in output data file (eg final.dat) and repopulates m with the parameters
 * @params *filename
 *  Filename containing parameters
 * @params *m
 *  Pointer to model to be populated
 * @return 1 on success, -1 on failure.
 */
int read_data_file(char *filename, struct Model * m) {
	/* This function reads a *.dx file as laid out in the docs.
	 * It reads in parameters to the struct m, and then calls
	 * read_relaxation and relevant other functions to import
	 * additional data.
	 *
	 * Function returns 1 if successful, -1 if failure.
	 */

	FILE * fp;
	char line[255];
	int len = 255;

	unsigned int k;
	for (k = 0; k < m->n_residues; k++)
		m->residues[k].ignore = 1;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}

	unsigned int params;
	switch (m->model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT: params= 5; break;
		case MOD_SMFT: params= 3; break;
		case MOD_DEMF: params = 4; break;
		case MOD_DEMFT: params= 6; break;
		case MOD_GAF: params = 8; break;
		case MOD_GAFT:params = 10; break;
		case MOD_EGAF: params = 6; break;
		case MOD_EGAFT: params = 8; break;
		default: params = 0; break;
	}
	if (m->or_variation == VARIANT_A)
		params += 3; // theta, phi

	int res;

	int length = 0;
	unsigned int length_d = 0;
	unsigned int i;


	while(fgets(line, len, fp)) {

		//3	0.780000	201.669690	1.461383e+01	9.604572e-01	0
		//int k = sscanf(line, "%d\t%lf\t%lf\t%n", &res, &S2NH, &chisq, &length);
		res--;
		m->residues[res].parameters = (Decimal *) malloc(sizeof(Decimal) * params);
		m->residues[res].ignore = 0;
		for (i = 0; i < params; i++) {
			k = sscanf(line + length, "%le\t%n", &(m->residues[res].parameters[i]), &length_d);
			//printf("%d: %lf\n", i, m->residues[res].parameters[i]);
			//m->residues[res].parameters[i] = (Decimal) param;
			//printf("K = %d\n", k);
			//printf("Param = %lf\n", parameters[i]);
			if (k <= 0) {
				break;
			}else {
				length += length_d;
			}
		}
	}

	/*int residue = 40;
	for (i = 0; i < params; i++) {
		printf("Parameter %d = %lf\n", i, m->residues[residue].parameters[i]);
	}*/

	fclose(fp);
	return 1;
}

/** 
 * outputs correlation function. Note that T and dT are both in nanoseconds (as this is the internal representation).
 * @params *fn
 *  Filename to write to
 * @params T
 *  Time (in ns) to run to.
 * @params dT
 *  Timestep (in ns) between outputs
 * @params taus
 *  Slow timescale in nanoseconds
 * @params S2s
 *  Slow motion order parameter
 * @params tauf
 *  Fast motion timescale in ns
 * @params S2f
 *  Fast motion order parameter
 * @return 1 on success, -1 on failure
 */
int write_correlation_function_emf(char * fn, Decimal T, Decimal dT, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f) {
	Decimal S2 = S2f * S2s;
	FILE * fp;

	fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("%s not found.\n", fn);
		return -1;
	}

	Decimal t;
	Decimal correl;
	int t_count;
	int t_count_max = (int) (T / dT);
	//for (t = 0; t < T; t += dT) {
	for (t_count = 0; t_count < t_count_max; t_count++) {
	    t = (Decimal) t_count * dT;
		correl = S2;
		correl += (1 - S2f) * expl(-t / tauf);
		correl += (S2f - S2) * expl(-t / taus);
		fprintf(fp, "%lf\t%lf\n", t, correl);
	}

	fclose(fp);
	return 1;
}

int write_spectral_density_emf(char *fn, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f) {
	FILE * fp;
	fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("%s not found.\n", fn);
		return -1;
	}
	
	Decimal w;
	Decimal v, slow, fast;
	Decimal w_start = 10 * T_DOWN;
	Decimal w_end = pow(10., 10) * T_DOWN;
	int n_steps = (int) (log(w_end / w_start) / log(2));
	int n;
	w = w_start;
	//for (w = pow(10, 1)*T_DOWN; w < pow(10., 10)*T_DOWN; w*=2) {
	for (n = 0; n < n_steps; n++) {
	    w = w * 2;
		/*v = (((1 - S2f) * tauf) / (1 + (w * tauf * w * tauf)));
		q = v + ((S2f * (1 - S2s) * taus) / (1 + (w * taus * w * taus)));
		
		l =  ( \
			((((1 - (Decimal) S2f)) * (Decimal) tauf)) \
			/ (1 + ((Decimal) w * (Decimal) tauf * (Decimal) w * (Decimal) tauf)) \
		);
		p = l +\
		(\
			(((Decimal) S2f) * (1 - (Decimal) S2s) * (Decimal) taus)\
			/ (1 + ((Decimal) w * (Decimal) taus * (Decimal) w * (Decimal) taus))\
		);*/
		v = J0(w, taus, S2s, tauf, S2f);
		//printf("%lf (%lf, %lf) (%lf, %lf) %lf\n", w, v, l, q, p, rt);
		fast = ( \
			((((1 - (Decimal) S2f)) * (Decimal) tauf)) \
			/ (1 + ((Decimal) w * (Decimal) tauf * (Decimal) w * (Decimal) tauf)) \
		);
		slow = (\
			(((Decimal) S2f) * (1 - (Decimal) S2s) * (Decimal) taus)\
			/ (1 + ((Decimal) w * (Decimal) taus * (Decimal) w * (Decimal) taus))\
		);
		v *= T_DOWN;
		slow *= T_DOWN;
		fast *= T_DOWN;
		printf("%lf\t%lf\n", w*T_UP, v);
		fprintf(fp, "%lf\t%le\t%le\t%le\n", w*T_UP, v, slow, fast);
	}
	fclose(fp);
	return 1;
}

int write_spectral_density_smf(char *fn, Decimal tau, Decimal S2) {
	FILE * fp;
	fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("%s not found.\n", fn);
		return -1;
	}
	
	Decimal w;
	Decimal v;
	Decimal w_start = pow(10, -6) * T_DOWN;
	Decimal w_end = pow(10, 6) * T_DOWN;
	int n_steps = (int) (log(w_end / w_start) / log(2));
	int n;
	w = w_start;
	for (n = 0; n < n_steps; n++) {
	//for (w = pow(10, -6)*T_DOWN; w < pow(10, 6)*T_DOWN; w*=2) {
	    w *= 2;
		v = J0(w, 0, 0, tau, S2) * T_DOWN;
		fprintf(fp, "%lf\t%le\n", w*T_UP, v);
	}
	fclose(fp);
	return 1;
}

/** 
 * outputs correlation function. Note that T and dT are both in nanoseconds (as this is the internal representation).
 * @params *fn
 *  Filename to write to
 * @params T
 *  Time (in ns) to run to.
 * @params dT
 *  Timestep (in ns) between outputs
 * @params tau
 *  Timescale in nanoseconds
 * @params S2
 *  Motion order parameter
 * @return 1 on success, -1 on failure
 */
int write_correlation_function_smf(char * fn, Decimal T, Decimal dT, Decimal tau, Decimal S2) {
	FILE * fp;

	fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("%s not found.\n", fn);
		return -1;
	}

	Decimal t;
	Decimal correl;
	int n, n_steps;
	n_steps = (int) (T / dT);
	for (n = 0; n < n_steps; n++) {
	//for (t = 0; t < T; t += dT) {
	    t = n * dT;
		correl = S2;
		correl += (1 - S2) * expl(-t / tau);
		fprintf(fp, "%lf\t%lf\n", t, correl);
	}

	fclose(fp);
	return 1;
}

/**
 * Runs correlation fitting. 
 * @opts -sSYSTEMFILE.dx
 *  Filename of system file (to obtain orientations, etc)
 * @opts -oOUTPUTFOLDER
 *  Folder to put all correlation plots into
 * @opts -fFINAL.DAT
 *  File containing final.dat filled params
 */
int main(int argc, char * argv[]) {

	/* Initialisation */
	start_time = time(0);
	srand((unsigned int)time(NULL));
	initialise_dwig(HALF_PI, Dwig);


	char system_file[255] = "";
	char data_file[255] = "";
	char output_folder[255] = "";
	unsigned int i;
	//int err_mod = 0;
	unsigned int bc_mod = 0, cor_mod = 0, sd_mod = 0;
	for (i = 1; i < (unsigned int) argc; i++) {
		//printf("%s\n", argv[i]);
		if (strncmp(argv[i], "-o", 2) == 0) {
			strcpy(output_folder, argv[i]+2);
			printf("OF: %s\n", output_folder);
		} else if (strncmp(argv[i], "-s", 2) == 0) {
			strcpy(system_file, argv[i]+2);
		} else if (strncmp(argv[i], "-f", 2) == 0) {
			strcpy(data_file, argv[i] + 2);
		} else if (strncmp(argv[i], "-bc", 3) == 0) {
			bc_mod = 1;
		} else if (strncmp(argv[i], "-co", 3) == 0) { 
			cor_mod = 1;
		} else if (strncmp(argv[i], "-de", 3) == 0) {
			sd_mod = 1;
		} else {
			printf("Please pass system file (-s), output folder (-o), and Dynamix results file (-f). Options -bc and -co.\n");
			return -1;
		}
	}

	if (strcmp(system_file, "") == 0 || strcmp(output_folder, "") == 0 || strcmp(data_file, "") == 0) {
		printf("Please pass system file (-s), output folder (-o), and Dynamix results file (-f). Options -bc and -co.\n");
		exit(-1);
	}



	struct Model m;
	int ret = read_system_file(system_file, &m);
	
	if (ret == -1) {
		printf("Reading system file failed.\n");
		free_all(&m);
		return -1;
	}

	ret = read_data_file(data_file, &m);
	printf("Ret: %d\n", ret);
	unsigned int params;
	switch (m.model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT: params= 5; break;
		case MOD_SMFT: params= 3; break;
		case MOD_DEMF: params = 4; break;
		case MOD_DEMFT: params= 6; break;
		case MOD_GAF: params = 8; break;
		case MOD_GAFT:params = 10; break;
		case MOD_EGAF: params = 6; break;
		case MOD_EGAFT: params = 8; break;
		default: params = 0; break;
	}
	if (m.or_variation == VARIANT_A)
		params += 3; // theta, phi


	
	int ignore = -1;
	#pragma omp parallel for 
	for (i = 0; i <m.n_residues; i++) {
//		Decimal Ea,Eas,Eaf, tauf,taus, tau, S2s, S2f;
		struct Residue * resid;
		Decimal * opts;
		unsigned int j;
//		Decimal S2NH, S2CN, S2CH, S2CC;
//		Decimal S2NHs, S2CNs, S2CHs, S2NHf, S2CHf, S2CNf, S2CCs, S2CCf;
//		Decimal temp;
		Decimal T = 200; // 1000
		Decimal dT = 0.01; // 0.01
		Decimal alpha, beta, gamma;
		char fn_NH[355];
		char fn_CN[355];
		char fn_CH[355];
		char fn_CC[355];
		char sd_NH[355];
		char sd_CN[355];
		char sd_CH[355];
		char sd_CC[355];
		char fn_BC[355];
		opts = (m.residues[i].parameters);
		printf("%d: %d\n", i, m.residues[i].ignore);
		if (m.residues[i].ignore == 1)
			continue;
		//printf("Hola. %d\n", sizeof(m.residues[i].parameters));
		for (j = 0; j < params; j++) {
			printf("%d %d %lf\n", i, j, m.residues[i].parameters[j]);
		}
		resid = &(m.residues[i]);
		if (cor_mod == 1 || sd_mod == 1) {

			sprintf(fn_NH, "%s/corr_%d_NH.dat", output_folder, i+1);
			sprintf(fn_CN, "%s/corr_%d_CN.dat", output_folder, i+1);
			sprintf(fn_CH, "%s/corr_%d_CH.dat", output_folder, i+1);
			sprintf(fn_CC, "%s/corr_%d_CCAp.dat", output_folder, i+1);
			sprintf(sd_NH, "%s/sd_%d_NH.dat", output_folder, i+1);
			sprintf(sd_CN, "%s/sd_%d_CN.dat", output_folder, i+1);
			sprintf(sd_CH, "%s/sd_%d_CH.dat", output_folder, i+1);
			sprintf(sd_CC, "%s/sd_%d_CCAp.dat", output_folder, i+1);
			printf("Hello\n");
			if (m.or_variation == VARIANT_A) {
				alpha = (Decimal) opts[params-3];
				beta = (Decimal) opts[params-2];
				gamma = (Decimal) opts[params-1];
				for (j = 0; j < N_OR; j++) {
					calculate_Y2(&(m.residues[i].orients[j]));
					rotate_Y2(&(m.residues[i].orients[j]), alpha, beta, gamma);
				}
			}

			int obts;
			struct BCParameters pars;
			obts = opts_to_bcpars(opts, &pars, m.model, &(m.residues[i]), &ignore);
            if (obts != 0)
                ERROR("Error converting opts to bcpars\n");

			Decimal taus = pars.taus, tauf = pars.tauf;
			if (pars.Eas != -1 && pars.Eaf != -1) {
			    taus = temp_tau(pars.taus, pars.Eas, 300);
			    tauf = temp_tau(pars.tauf, pars.Eaf, 300);
			}
			if (cor_mod == 1) {
				if (m.model == MOD_SMF || m.model == MOD_SMFT) {
					write_correlation_function_smf(fn_NH, T, dT, tauf, pars.S2NHs * pars.S2NHf);
					write_correlation_function_smf(fn_CH, T, dT, tauf, pars.S2CHs * pars.S2CHf);
					write_correlation_function_smf(fn_CN, T, dT, tauf, pars.S2CNs * pars.S2CNf);
					write_correlation_function_smf(fn_CC, T, dT, tauf, pars.S2CCs * pars.S2CCf);
				} else {
					write_correlation_function_emf(fn_NH, T, dT, taus, pars.S2NHs, tauf, pars.S2NHf);
					write_correlation_function_emf(fn_CH, T, dT, taus, pars.S2CHs, tauf, pars.S2CHf);
					write_correlation_function_emf(fn_CN, T, dT, taus, pars.S2CNs, tauf, pars.S2CNf);
					write_correlation_function_emf(fn_CC, T, dT, taus, pars.S2CCs, tauf, pars.S2CCf);
				}
			}
			if (sd_mod == 1) {
				if (m.model == MOD_SMF || m.model == MOD_SMFT) {
					write_spectral_density_smf(sd_NH, tauf, pars.S2NHs * pars.S2NHf);
					write_spectral_density_smf(sd_CH, tauf, pars.S2CHs * pars.S2CHf);
					write_spectral_density_smf(sd_CN, tauf, pars.S2CNs * pars.S2CNf);
					write_spectral_density_smf(sd_CC, tauf, pars.S2CCs * pars.S2CCf);
				} else {
					write_spectral_density_emf(sd_NH, taus, pars.S2NHs, tauf, pars.S2NHf);
					write_spectral_density_emf(sd_CH, taus, pars.S2CHs, tauf, pars.S2CHf);
					write_spectral_density_emf(sd_CN, taus, pars.S2CNs, tauf, pars.S2CNf);
					write_spectral_density_emf(sd_CC, taus, pars.S2CCs, tauf, pars.S2CCf);
				}
			}
		}
		if (bc_mod == 1) {
			sprintf(fn_BC, "%s/backcalc_%d.dat", output_folder, i+1);
			back_calculate(opts, resid, &m, fn_BC, params);


			// int back_calculate(Decimal * opts, struct Residue * resid, unsigned int model, unsigned int or_variations, char *filename, unsigned int params) {

		}
	}
	//int write_correlation_function_smf(char * fn, Decimal T, Decimal dT, Decimal tau, Decimal S2) {


	//fclose(fp);
	free_all(&m);
	return 1;
}
