/**
 * @file read_data.c
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/**
 * Reads residue specific data from datafile. \n
 * Data file being read should be set out in columns, with spaces between them.\n
 * \t Column one should contain the residue number as an integer, indexed from 1.\n
 * \t Column two should contain the data of note.\n
 * \t Column three should contain the error (or 0).\n
 * @note Comments may be inserted by preceeding a line with '%'\n
 * @param m
 *  Pointer to model
 * @param filename
 *  File to read
 * @param dt
 *  Data type (one of DATA_S2, DATA_CSISOC, DATA_CSISON)
 * @return 1 if successful, -1 else.
 */
int read_resid_data(struct Model *m, char *filename, int dt) {
	FILE *fp;
	char line[255];
	size_t len=255;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	int resid;
	float val;
	float err;
	
	while(fgets(line, len, fp)) {
		if (line[0] == '%')
			continue; // comment
		int k = sscanf(line, "%d %f %f\n", &resid, &val, &err);
		if (k != 3)
			printf("Error reading '%s'\n", line);
		resid = resid - 1; // 0 indexed in C
		switch (dt) {
			case DATA_S2:
				m->residues[resid].S2_dipolar = val;
				if (val == -1) m->residues[resid].ignore = 1;
				break;
			case DATA_CSISON:
				m->residues[resid].csisoN = val;
				/* 15N CSA parametrization according to isotropic chemical shift
				 * Wylie, B. J. et al.; Rienstra, C. M. Proc. Natl. Acad. Sci. U. S. A. 2011, 108, 16974.%from Rienstra, PNAS 2011.:
				 */
				m->residues[resid].csaN[0] = 1.1283 * val + 93.77;
				m->residues[resid].csaN[1] = 1.0086 * val - 42.475;
				m->residues[resid].csaN[2] = 0.8631 * val - 51.295;
				break;
			case DATA_CSISOC:
				m->residues[resid].csisoC = val;
				/* 13C' CSA parametrization according to isotropic chemical shift from
				 * Wylie, B. J. et al.; Rienstra, C. M. J. Am. Chem. Soc. 2007, 129, 5318.*/
				m->residues[resid].csaC[0] = 0.24 * val + 200;
				m->residues[resid].csaC[1] = 2.82 * val - 305;
				m->residues[resid].csaC[2] = 96.5;
				break;
			default:break;
		}
	}
	fclose(fp);
	return 1;
}

/**
 * Reads peptide plane orientation data. Should be in a file with 3 spaced columns;\n
 * \t Column 1 should contain the residue number indexed from 1\n
 * \t Column 2 should contain the theta angle (in radians)\n
 * \t Column 3 should contain the phi angle (in radians)\n
 * @note Comments may be inserted by preceeding a line with '%'\n
 * @param m
 *  Pointer to model
 * @param filename
 *  Filename
 * @param orient
 *  One of OR_* denoting which orientations are being read.
 * @return 1 if successful, else -1.
 */
int read_pp(struct Model *m, char *filename, int orient) {
	FILE *fp;
	char line[255];
	size_t len=255;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	int resid;
	float theta;
	float phi;
	while(fgets(line, len, fp)) {
		if (line[0] == '%')
			continue; // comment
		int k = sscanf(line, "%d %f %f\n", &resid, &theta, &phi);
		if (k != 3)
			printf("Error reading line '%s'\n", line);
		resid = resid -1; // index from 0
		m->residues[resid].orients[orient].theta = theta;
		m->residues[resid].orients[orient].phi = phi;
		calculate_Y2(&(m->residues[resid].orients[orient]));
	}
	return 1;
}

/**
 * Reads relaxation data. Relaxation data should be set out with two components.\n
 * The file should begin with the keys 'FIELD', 'WR', 'W1', 'TYPE' denoting the form of the data, eg\n
 * \t FIELD = 700\n
 * \t WR = 50000\n
 * \t W1 = 0\n
 * \t TEMP = 300\n
 * \t TYPE = 15NR1\n
 * Each parameter will default to -1 so if one is not set this may cause errors (though a warning will be given). \n
 * TYPE may be one of [15NR1, 15NR1p, 13CR1, 13CR1p] currently, though it should be noted that 13C has not been tested.\n
 * @warning There must be a space on each size of the equals sign (eg, TEMP = 300 is good, TEMP=300 is bad)\n
 * Once this head section is complete, it should be followed by "#DATA".\n
 * After this, data should be split into three columns;\n
 * \t Column 1: residue number (starting from 1)\n
 * \t Column 2: Relaxation rate (in s-1)\n
 * \t Column 3: Relaxation error, two standard deviations (in s-1)\n
 * @note Comments may be inserted by preceeding a line with '%'\n
 * @param m
 *  Pointer to model
 * @param filename
 *  Filename
 * @return 1 if successful, else -1.
 */
int read_relaxation_data(struct Model *m, char *filename) {
	FILE *fp;
	char line[255];
	size_t len=255;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	int mode = 0;
	char key[255], val[255];
	int i;
	float field = -1; // in MHz
	float wr = -1; // in Hz
	float w1 = -1; // in Hz
	int type = -1;
	float T = -1; // in Kelvin
	int resid;
	int rel = -1;
	float R, Re;
	while(fgets(line, len, fp)) {
		if (line[0] == '%')
			continue; // comment
		if (strcmp(line, "#DATA\n") == 0) {
			for (i = 0; i < m->n_residues; i++) {
				/* first check we have enough memory */
				if (m->residues[i].lim_relaxation - m->residues[i].n_relaxation < 5) {
					printf("Oof\n");
					m->residues[i].lim_relaxation += N_RELAXATION;
					m->residues[i].relaxation = (struct Relaxation *) realloc(m->residues[i].relaxation, sizeof(struct Relaxation) * (m->residues[i].lim_relaxation));
					if (m->residues[i].relaxation == NULL) {
						printf("Pointer loss\n");
						fclose(fp); // clean up on the way out
						exit(-1);
					}
				}
				rel = m->residues[i].n_relaxation;
				m->residues[i].n_relaxation++;
				if (field == -1)
					printf("Note: Field = 0 in %s\n", filename);
				if (wr == -1)
					printf("Note: wr = 0 in %s\n", filename);
				if (w1 == -1)
					printf("Note: w1 = 0 in %s\n", filename);
				if (type == -1)
					printf("Note: type = 0 in %s\n", filename);
				if (T == -1)
					printf("Note: T = 0 in %s\n", filename);

				m->residues[i].relaxation[rel].field = field;
				m->residues[i].relaxation[rel].wr = wr;
				m->residues[i].relaxation[rel].w1 = w1;
				m->residues[i].relaxation[rel].type = type;
				m->residues[i].relaxation[rel].T = T;

			}
			mode = 1;
			continue;
		}
		if (mode == 0) {
			sscanf(line, "%255s = %255s\n", key, val);
			if (strcmp(key, "FIELD") == 0)
				field = atof(val);
			else if (strcmp(key, "WR") == 0)
				wr = atof(val);
			else if (strcmp(key, "W1") == 0)
				w1 = atof(val);
			else if (strcmp(key, "TEMP") == 0)
				T = atof(val);
			else if (strcmp(key, "TYPE") == 0) {
				if (strcmp(val, "15NR1") == 0)
					type = R_15NR1;
				else if (strcmp(val, "15NR1p") == 0)
					type = R_15NR1p;
				else if (strcmp(val, "13CR1") == 0)
					type = R_13CR1;
				else if (strcmp(val, "13CR1p") == 0)
					type = R_13CR1p;
			} else {
				printf("Parameter %s unknown.\n", key);
				fclose(fp);
				free_all(m);
				return -1;
			}
		} else if (mode == 1) {
			int k = sscanf(line, "%d %f %f\n", &resid, &R, &Re);
			if (k != 3) {
				printf("Read error.\n");
			}/* else {
				printf("%s : %d, %f %f\n", line, resid, R, Re);
			}*/
			//printf("%d, %d, %f, %f\n", rel, resid, R, Re);
			// index from 0
			/*if (R == -1) {
				//printf("%s %d\n", line, m->residues[resid-1].n_relaxation);
				m->residues[resid-1].n_relaxation = m->residues[resid-1].n_relaxation - 1;
				if (m->residues[resid-1].n_relaxation < 0)
					m->residues[resid-1].n_relaxation = 0;
				continue;
			} else {*/
			m->residues[resid-1].relaxation[rel].R = R;
			m->residues[resid-1].relaxation[rel].Rerror = Re; // note indexing from 0
			//}
			
		}
	}
	fclose(fp);
	return 1;
}

/**
 * Reads system file. System file should begin with key value pairs laid out as\n
 * \t KEY = VALUE\n
 * With spaces either side of the equals. Keys which may be used are as follows;\n
 * \t MODEL (one of SMF, SMFT, EMF, EMFT, GAF, GAFT)\n
 * \t S2DIP (file containing dipolar order parameters)\n
 * \t CSISON (file containing isotropic chemical shifts for N)\n
 * \t CSISOC (file containign isotropic chemical shifts for C)\n
 * \t N_RESIDUES (number of residues in protein)\n
 * \t OUTPUT (directory for output data)\n
 * \t N_ITER (number of iterations for main fitting) \n
 * \t IGNORE (may be repeated, once for each residue to ignore.)\n
 * \t N_ERROR_ITER (number of iterations for error calculation) \n
 * \t OR_* (file containing orientation data. May be one of OR_NH, OR_NC, OR_NCA, OR_NCSAxx, OR_NCSAyy, OR_NCSAzz, OR_CCAp, OR_CCAc, OR_CN, OR_CNH, OR_CCSAxx, OR_CCAyy, OR_CCAzz. Should be set out as in read_pp().)\n
 * \t NTHREADS (number of parallel threads to spawn)\n
 * This should then be followed by "#RELAXATION". After this, the files containing relaxation data (as set out in read_relaxation_data()) should be listed.
 * 
 * @note Comments may be inserted by preceeding a line with '%'
 * @warning There must be a space on each size of the equals sign (eg, TEMP = 300 is good, TEMP=300 is bad)\n
 * @param filename
 *  filename
 * @param m
 *  Pointer to model
 * @return 1 if successful, else -1.
 */
int read_system_file(char *filename, struct Model * m) {
	/* This function reads a *.dx file as laid out in the docs.
	 * It reads in parameters to the struct m, and then calls
	 * read_relaxation and relevant other functions to import
	 * additional data.
	 *
	 * Function returns 1 if successful, -1 if failure.
	 */

	FILE * fp;
	char line[255];
	size_t len = 255;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	int mode = 0;
	char key[255];
	char val[255];
	char s2d[255] = "";
	char csisoN[255] = "";
	char csisoC[255] = "";
	int n_resid = -1;
	char pp_orient[14][255];
	int i;
	for (i = 0; i < 14; i++) {
		strcpy(pp_orient[i], "");
	}
	strcpy(m->outputdir, "./");
	m->n_iter = 1;
	int to_ignore[50];
	int ig = 0;

	/* Initialise */
	m->max_func_evals = 20000;
	m->max_iter = 20000;
	m->n_iter = -1;
	m->n_error_iter = -1;
	m->model = -1;
	m->n_residues = -1;
	m->nthreads = 1;


	while(fgets(line, len, fp)) {
		if (line[0] == '%')
			continue; // comment
		if (strcmp(line, "#RELAXATION\n") == 0) {
			if (n_resid == -1) {
				printf("Error: Number of residues must be defined\n");
				return -1;
			}

			m->residues = (struct Residue *) malloc(sizeof(struct Residue) * n_resid);
			int i;
			for (i = 0; i < n_resid; i++) {
				m->residues[i].relaxation = (struct Relaxation *) malloc(sizeof(struct Relaxation) * N_RELAXATION);
				m->residues[i].lim_relaxation = N_RELAXATION;
				m->residues[i].n_relaxation = 0;
				m->residues[i].ignore = 0;
			}

			for (i = 0; i < ig; i++) {
				m->residues[to_ignore[i] - 1].ignore = 1;
			}

			int t=0;
			if (strcmp(s2d, "") != 0)
				t += read_resid_data(m, s2d, DATA_S2);
			else
				t++;
			if (strcmp(csisoN, "") != 0)
				t += read_resid_data(m, csisoN, DATA_CSISON);
			else
				t++;
			if (strcmp(csisoC, "") != 0)
				t += read_resid_data(m, csisoC, DATA_CSISOC);
			else
				t++;

			if (t != 3) {
				printf("Error: Error in reading one of S2, csisoN, csisoC\n");
				return -1;
			}
			t = 0;
			for (i = 0; i < 14; i++) {
				if (strcmp(pp_orient[i], "") != 0)
					t += read_pp(m, pp_orient[i], i);
				else
					t++;
			}
			if (t != 14) {
				printf("Error: Error in reading orientations\n");
				return -1;
			}
			mode = 1;
			continue;
		}
		if (mode == 0) {
			if (line[0] == '#')
				continue;
			sscanf(line, "%255s = %255s\n", key, val);
			if (strcmp(key, "MODEL") == 0) {
				if (strcmp(val, "SMF") == 0)
					m->model = MOD_SMF;
				else if (strcmp(val, "EMF") == 0)
					m->model = MOD_EMF;
				else if (strcmp(val, "EMFT") == 0)
					m->model = MOD_EMFT;
				else if (strcmp(val, "SMFT") == 0)
					m->model = MOD_SMFT;
				else if (strcmp(val, "GAF") == 0)
					m->model = MOD_GAF;
				else if (strcmp(val, "GAFT") == 0)
					m->model = MOD_GAFT;
				else if (strcmp(val, "DEMF") == 0)
					m->model = MOD_DEMF;
				else if (strcmp(val, "DEMFT") == 0)
					m->model = MOD_DEMFT;
				else {
					printf("Model %s unknown.\n", val);
					return -1;
				}
			} else if (strcmp(key, "S2DIP") == 0) {
				strcpy(s2d, val);
			} else if (strcmp(key, "CSISON") == 0) {
				strcpy(csisoN, val);
			} else if (strcmp(key, "CSISOC") == 0) {
				strcpy(csisoC, val);
			} else if (strcmp(key, "MAXFUNCEVALS") == 0) {
				m->max_func_evals = atoi(val);
			} else if (strcmp(key, "MAXITER") == 0) {
				m->max_iter = atoi(val);
			} else if (strcmp(key, "N_RESIDUES") == 0){
				m->n_residues = atoi(val);
				n_resid = atoi(val);
			} else if (strcmp(key, "N_ITER") == 0) {
				m->n_iter = atoi(val);
			} else if (strcmp(key, "OUTPUT") == 0) {
				strcpy(m->outputdir, val);
			} else if (strcmp(key, "IGNORE") == 0) {
				to_ignore[ig] = atoi(val);
				ig++;
			} else if (strcmp(key, "N_ERROR_ITER") == 0){
				m->n_error_iter = atoi(val);
			} else if (strcmp(key, "OR_NH") == 0) {
				strcpy(pp_orient[OR_NH], val);
			} else if (strcmp(key, "OR_NC") == 0) {
				strcpy(pp_orient[OR_NC], val);
			} else if (strcmp(key, "OR_NCA") == 0) {
				strcpy(pp_orient[OR_NCA], val);
			} else if (strcmp(key, "OR_NCSAxx") == 0) {
				strcpy(pp_orient[OR_NCSAxx], val);
			} else if (strcmp(key, "OR_NCSAyy") == 0) {
				strcpy(pp_orient[OR_NCSAyy], val);
			} else if (strcmp(key, "OR_NCSAzz") == 0) {
				strcpy(pp_orient[OR_NCSAzz], val);
			} else if (strcmp(key, "OR_CCAp") == 0) {
				strcpy(pp_orient[OR_CCAp], val);
			} else if (strcmp(key, "OR_CCAc") == 0) {
				strcpy(pp_orient[OR_CCAc], val);
			} else if (strcmp(key, "OR_CN") == 0) {
				strcpy(pp_orient[OR_CN], val);
			} else if (strcmp(key, "OR_CNH") == 0) {
				strcpy(pp_orient[OR_CNH], val);
			} else if (strcmp(key, "OR_CH") == 0) {
				strcpy(pp_orient[OR_CH], val);
			} else if (strcmp(key, "OR_CCSAxx") == 0) {
				strcpy(pp_orient[OR_CCSAxx], val);
			} else if (strcmp(key, "OR_CCSAyy") == 0) {
				strcpy(pp_orient[OR_CCSAyy], val);
			} else if (strcmp(key, "OR_CCSAzz") == 0) {
				strcpy(pp_orient[OR_CCSAzz], val);
			} else if (strcmp(key, "NTHREADS") == 0) {
				m->nthreads = atoi(val);
				if (m->nthreads <= 0) {
					printf("Number of threads must be greater than 0.\n");
					exit(-1);
				}
			}

			//printf("%s: %s\n", key, val);
		} else if (mode == 1) {
			sscanf(line, "%255s\n", key);
			//printf("Read file: %s\n", key);
			//int read_relaxation_data(struct Model *m, char *filename) {

			read_relaxation_data(m, key);
		}
	}
	
	for (i = 0; i < m->n_residues; i++) {
		if (m->residues[i].n_relaxation < 20)
			m->residues[i].ignore = 1;
	}
	fclose(fp);

	/*Check requirements */
	if (m->n_iter < 0) {
		printf("Please provide N_ITER\n");
		return -1;
	}
	if (m->model < 0) {
		printf("Please provide MODEL\n");
		return -1;
	}
	if (m->n_residues < 0){
		printf("Please provide residue count\n");
		return -1;
	}

	return 1;


	/*	int max_func_evals, max_iter;
	int model;
	struct Residue * residues;
	int n_residues;*/

}

int print_system(struct Model *m, char *filename) {
	FILE * fp;
	fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	//fprintf(fp, "MFE: %d\nMI:  %d\n", m->max_func_evals, m->max_iter);
	fprintf(fp, "Model: %d\nN_Residues: %d\n", m->model, m->n_residues);
	int i, j;
	struct Relaxation *r;
	for (i = 0; i < m->n_residues; i++) {
		fprintf(fp, "=== Residue %d ===\n", i+1);
		if (m->residues[i].ignore == 1)
			fprintf(fp, "--- IGNORING ---\n");

		fprintf(fp, "\tS2 = %f\n\tCSISON = %f\n\tCSISOC = %f\n", m->residues[i].S2_dipolar, m->residues[i].csisoN, m->residues[i].csisoC);
		fprintf(fp, "\tCSAN: [%f, %f, %f]\n", m->residues[i].csaN[0], m->residues[i].csaN[1], m->residues[i].csaN[2]);
		fprintf(fp, "\tCSAC: [%f, %f, %f]\n", m->residues[i].csaC[0], m->residues[i].csaC[1], m->residues[i].csaC[2]);
		fprintf(fp, "\tOrients;\n");
		for (j = 0; j < 14; j++) {
			//if (m->residues[i].orients[j] != NULL)
			fprintf(fp, "\t\t%d: %f, %f\n", j, m->residues[i].orients[j].theta, m->residues[i].orients[j].phi);
		}
		fprintf(fp, "\tRelaxation Constraints: %d\n", m->residues[i].n_relaxation);
		for (j = 0; j < m->residues[i].n_relaxation; j++) {
			r = &(m->residues[i].relaxation[j]);
			fprintf(fp, "\t\t%d: %d [%f, %f, %f, %f] %f +- %f\n", j, r->type, r->field, r->wr, r->w1, r->T, r->R, r->Rerror);
		}
	}
	fclose(fp);
	return 1;

}