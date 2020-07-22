#include <stdlib.h>
#define MOD_SMF 		0
#define MOD_EMF 		1
#define MOD_EMFT		2
#define MOD_SMFT 		3
#define MOD_GAF 		4
#define MOD_GAFT 		5

#define DATA_S2			0
#define DATA_CSISON		1
#define DATA_CSISOC		2

#define R_15NR1 		0
#define R_15NR1p		1
#define R_13CR1 		2
#define R_13CR1p		3

#define OR_NH			0
#define OR_NC			1
#define OR_NCA			2
#define OR_NCSAxx		3
#define OR_NCSAyy 		4
#define OR_NCSAzz		5
#define OR_NCCAp		6
#define OR_NCCAc		7
#define OR_CN			8
#define OR_CNH			9
#define OR_CH			10
#define OR_CCSAxx		11
#define OR_CCSAyy		12
#define OR_CCSAzz		13

#define THETA 			0
#define PHI				1

#define N_RELAXATION	50
#define NTHREADS		4	
#define MIN_VAL			10000000

#define RYD 			8.3144621

/* These constants are extracted from NHCO_3GAFSquaredEa, using the bond lengths defined in there */
#define D_NH			72038.41107
#define D_CC			13467.66052
#define D_CH			22362.47724
#define D_CHr			31491.71046
#define D_CN			8175.2
#define D_CaN			6180.1
#define D_HNr			13108.32273

#define MODE_15N		0
#define MODE_13C		1

struct Model {
	int max_func_evals, max_iter, n_iter, n_error_iter;
	char outputdir[255];
	int model;
	struct Residue * residues;
	int n_residues;
};

struct Residue {
	float S2_dipolar;
	float csisoN;
	float csisoC;
	float csaN[3];
	float csaC[3];
	float orients[14][2]; // theta,phi
	struct Relaxation * relaxation;
	struct Relaxation * temp_relaxation;
	long double ** error_params;
	
	int n_relaxation;
	int lim_relaxation;
	int ignore;
	double min_val; 
	long double * parameters;
	long double * errors_std;
	long double * errors_mean;
};

struct Relaxation {
	float field; // in MHz
	float wr; // in Hz
	float w1; // in Hz
	int type;
	float R;
	float Rerror;
	float T; // in Kelvin
};

struct rrargs {
	int i;
	struct Residue * resid;
	int model;
	int n_iter;
	char outputdir[255];
};

void free_all(struct Model *m) {
	int res, k;
	int params;
	switch (m->model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT:params = 5; break;
		case MOD_SMFT:params = 3; break;
		default: params = 0; break;
	}
	
	for (res = 0; res < m->n_residues; res++) {
		free(m->residues[res].relaxation);
		free(m->residues[res].parameters);
		if (m->residues[res].error_params != NULL) {
			for (k = 0; k < params; k++) {
				if (m->residues[res].error_params[k] != NULL)
					free(m->residues[res].error_params[k]);
			}
		}
		if (m->residues[res].errors_mean != NULL)
			free(m->residues[res].errors_mean);
		if (m->residues[res].errors_std != NULL)
			free(m->residues[res].errors_std);
	}
	free(m->residues);
	
	return;
}
