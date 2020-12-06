/**
 * @file correlation.c
 *  Calculates correlation functions, compiled separately.
 */

#include <stdio.h>
#include "datatypes.h"
#include "read_data.h"

#include "models/smf.h"
#include "models/emf.h"
#include "models/gaf.h"
#include "models/egaf.h"
#include "chisq.h"
#include "crosen.h" // implementation of Nelder-Mead simplex algorithm
#include "errors.h"
#include <stdlib.h>  
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

	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}

	unsigned int params = 0;
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
	float S2NH, chisq;
	int length = 0;
	int length_d = 0;
	unsigned int i;


	while(fgets(line, len, fp)) {

		//3	0.780000	201.669690	1.461383e+01	9.604572e-01	0
		int k = sscanf(line, "%d\t%f\t%f\t%n", &res, &S2NH, &chisq, &length);
		res--;
		m->residues[res].parameters = (long double *) malloc(sizeof(long double) * params);

		for (i = 0; i < params; i++) {
			k = sscanf(line + length, "%Le\t%n", &(m->residues[res].parameters[i]), &length_d);
			//m->residues[res].parameters[i] = (long double) param;
			//printf("K = %d\n", k);
			//printf("Param = %f\n", parameters[i]);
			if (k <= 0)
				break;
			else
				length += length_d;
		}
	}

	/*int residue = 40;
	for (i = 0; i < params; i++) {
		printf("Parameter %d = %Lf\n", i, m->residues[residue].parameters[i]);
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
int write_correlation_function_emf(char * fn, double T, double dT, long double taus, long double S2s, long double tauf, long double S2f) {
	long double S2 = S2f * S2s;
	FILE * fp;

	fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("%s not found.\n", fn);
		return -1;
	}

	double t;
	double correl = 0;
	for (t = 0; t < T; t += dT) {
		correl = S2;
		correl += (1 - S2f) * expl(-t / tauf);
		correl += (S2f - S2) * expl(-t / taus);
		fprintf(fp, "%lf\t%lf\n", t, correl);
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
int write_correlation_function_smf(char * fn, double T, double dT, long double tau, long double S2) {
	FILE * fp;

	fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("%s not found.\n", fn);
		return -1;
	}

	double t;
	double correl = 0;
	for (t = 0; t < T; t += dT) {
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
	unsigned int bc_mod = 0, cor_mod = 0;
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
	long double * opts;
	unsigned int j;
	long double S2NH=0, S2CN=0, S2CH=0;
	double S2NHs, S2CNs, S2CHs, S2NHf, S2CHf, S2CNf;
	long double temp = 300;
	double T = 200; // 1000
	double dT = 0.01; // 0.01
	double alpha, beta, gamma;
	char fn_NH[355];
	char fn_CN[355];
	char fn_CH[355];
	char fn_BC[355];

	unsigned int params = 0;
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


	long double Ea=0,Eas=0,Eaf=0, tauf=0,taus=0, tau=0, S2s=0, S2f=0;;
	struct Residue * resid;
	for (i = 0; i < m.n_residues; i++) {
		opts = (m.residues[i].parameters);
		resid = &(m.residues[i]);
		if (cor_mod == 1) {
			
			sprintf(fn_NH, "%s/corr_%d_NH.dat", output_folder, i+1);
			sprintf(fn_CN, "%s/corr_%d_CN.dat", output_folder, i+1);
			sprintf(fn_CH, "%s/corr_%d_CH.dat", output_folder, i+1);

			if (m.or_variation == VARIANT_A) {
				alpha = (double) opts[params-3];
				beta = (double) opts[params-2];
				gamma = (double) opts[params-1];
				for (j = 0; j < N_OR; j++) {
					calculate_Y2(&(m.residues[i].orients[j]));
					rotate_Y2(&(m.residues[i].orients[j]), alpha, beta, gamma);
				}
			}

			if (m.model == MOD_SMF || m.model == MOD_SMFT) {
				tau = opts[0];
				S2NH = (double) opts[1];
				S2CN = (double) opts[1];
				S2CH = (double) opts[1];
				if (m.model == MOD_SMFT) {
					Ea = opts[2];
					tau *= expl(Ea / (RYD * temp));
				}
				if (S2NH < 0 || tau < 0)
					continue;
				
			} else if (m.model == MOD_EMF || m.model == MOD_EMFT || m.model == MOD_DEMF || m.model == MOD_DEMFT) {
				taus = opts[0];
				S2s = opts[1];
				tauf = opts[2];
				S2f = m.residues[i].S2NH / S2s;
				if (m.model == MOD_EMFT) {
					Eas = opts[3];
					Eaf = opts[4];
				}  else if (m.model == MOD_DEMFT) {
					Eas = opts[4];
					Eaf = opts[5];
				}
				taus *= expl(Eas / (RYD * temp));
				tauf *= expl(Eaf / (RYD * temp));
				if (m.model == MOD_DEMF || m.model == MOD_DEMFT) {
					S2f = opts[3];
				}
				S2NHs = (double) S2s;
				S2CHs = (double) S2s;
				S2CNs = (double) S2s;
				S2NHf = (double) S2f;
				S2CHf = (double) S2f;
				S2CNf = (double) S2f;
			} else if (m.model == MOD_GAF || m.model == MOD_GAFT) {
				// need to perform reorientation before.
				taus = opts[0];
				tauf = opts[1];
				long double sigs[3] = {opts[2], opts[3], opts[4]};
				long double sigf[3] = {opts[5], opts[6], opts[7]};
				if (m.model == MOD_GAFT) {
					Eas = opts[8];
					Eaf = opts[9];
					taus *= expl(Eas / (RYD * temp));
					tauf *= expl(Eaf / (RYD * temp));
				}
				struct Orient *Os[] = {&(resid->orients[OR_NH]), &(resid->orients[OR_CNH]), &(resid->orients[OR_CN])};
				double *S2s[] = {&S2NHs, &S2CHs, &S2CNs};
				double *S2f[] = {&S2NHf, &S2CHf, &S2CNf};
				GAF_S2(sigs, Os, Os, S2s, 3, MODE_REAL);
				GAF_S2(sigf, Os, Os, S2f, 3, MODE_REAL);
			} else if (m.model == MOD_EGAF || m.model == MOD_EGAFT) {
				// need to perform reorientation before.
				taus = opts[0];
				tauf = opts[1];
				long double sigs[3] = {opts[2], opts[3], opts[4]};
				S2f = opts[5];
				if (m.model == MOD_GAFT) {
					Eas = opts[6];
					Eaf = opts[7];
					taus *= expl(Eas / (RYD * temp));
					tauf *= expl(Eaf / (RYD * temp));
				}
				struct Orient *Os[] = {&(resid->orients[OR_NH]), &(resid->orients[OR_CNH]), &(resid->orients[OR_CN])};
				double *S2s[] = {&S2NHs, &S2CHs, &S2CNs};
				GAF_S2(sigs, Os, Os, S2s, 3, MODE_REAL);
				S2NHf = (double) S2f;
				S2CNf = (double) S2f;
				S2CHf = (double) S2f;
			}

			if (m.model == MOD_SMF || m.model == MOD_SMFT) {
				write_correlation_function_smf(fn_NH, T, dT, tau, S2NH);
				write_correlation_function_smf(fn_CH, T, dT, tau, S2CH);
				write_correlation_function_smf(fn_CN, T, dT, tau, S2CN);
			} else {
				write_correlation_function_emf(fn_NH, T, dT, taus, S2NHs, tauf, S2NHf);
				write_correlation_function_emf(fn_CH, T, dT, taus, S2CHs, tauf, S2CHf);
				write_correlation_function_emf(fn_CN, T, dT, taus, S2CNs, tauf, S2CNf);
			}
		}
		if (bc_mod == 1) {
			sprintf(fn_BC, "%s/bc_%d.csv", output_folder, i+1);
			back_calculate(opts, resid, &m, fn_BC, params);


			// int back_calculate(long double * opts, struct Residue * resid, unsigned int model, unsigned int or_variations, char *filename, unsigned int params) {

		}
	}
	//int write_correlation_function_smf(char * fn, double T, double dT, long double tau, long double S2) {


	//fclose(fp);
	free_all(&m);
	return 1;
}
