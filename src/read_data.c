#include <stdio.h>
#include <string.h>
#include <stdlib.h>


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
	int k;
	while(fgets(line, len, fp)) {
		if (line[0] == '%')
			continue; // comment
		k = sscanf(line, "%d %f %f\n", &resid, &val, &err);
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
	int k;
	while(fgets(line, len, fp)) {
		if (line[0] == '%')
			continue; // comment
		k = sscanf(line, "%d %f %f\n", &resid, &theta, &phi);
		resid = resid -1; // index from 0
		m->residues[resid].orients[orient][THETA] = theta;
		m->residues[resid].orients[orient][PHI] = phi;
	}
	return 1;
}

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
			} else {
				printf("Parameter %s unknown.\n", key);
				fclose(fp);
				free_all(m);
				return -1;
			}
		} else if (mode == 1) {
			sscanf(line, "%d %f %f\n", &resid, &R, &Re);
			//printf("%d, %d, %f, %f\n", rel, resid, R, Re);
			// index from 0
			m->residues[resid-1].relaxation[rel].R = R;
			if (R == -1)
				m->residues[resid-1].ignore = 1;
			m->residues[resid-1].relaxation[rel].Rerror = Re; // note indexing from 0
		}
	}
	fclose(fp);
	return 1;
}

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
	for (i = 0; i < 10; i++) {
		strcpy(pp_orient[i], "");
	}
	strcpy(m->outputdir, "./");
	m->n_iter = 1;
	int to_ignore[50];
	int ig = 0;
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
			} else if (strcmp(key, "OR_NCCAp") == 0) {
				strcpy(pp_orient[OR_NCCAp], val);
			} else if (strcmp(key, "OR_NCCAc") == 0) {
				strcpy(pp_orient[OR_NCCAc], val);
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
			}

			//printf("%s: %s\n", key, val);
		} else if (mode == 1) {
			sscanf(line, "%255s\n", key);
			//printf("Read file: %s\n", key);
			//int read_relaxation_data(struct Model *m, char *filename) {

			read_relaxation_data(m, key);
		}
	}
	fclose(fp);
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
	fprintf(fp, "MFE: %d\nMI:  %d\n", m->max_func_evals, m->max_iter);
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
			fprintf(fp, "\t\t%d: %f, %f\n", j, m->residues[i].orients[j][0], m->residues[i].orients[j][1]);
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

