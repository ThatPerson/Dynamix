#include <stdlib.h>
#define MOD_SMF 	0
#define MOD_EMF 	1
#define MOD_EMFT	2
#define MOD_SMFT 	3
#define MOD_GAF 	4
#define MOD_GAFT 	5

#define DATA_S2		0
#define DATA_CSISON	1
#define DATA_CSISOC	2

#define R_15NR1 	0
#define R_15NR1p	1

#define OR_NH		0
#define OR_NC		1
#define OR_NCA		2
#define OR_NCSA		3
#define OR_NCCAp	4
#define OR_NCCAc	5
#define OR_CN		6
#define OR_CNH		7
#define OR_CH		8

#define THETA 		0
#define PHI			1

#define N_RELAXATION 50

struct Model {
	int max_func_evals, max_iter, n_iter;
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
	float orients[10][2]; // theta,phi
	struct Relaxation * relaxation;
	int n_relaxation;
	int lim_relaxation;
	int ignore;
	double min_val; 
	long double parameters[20];
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
	int res;
	for (res = 0; res < m->n_residues; res++) {
		free(m->residues[res].relaxation);
	}
	free(m->residues);
	return;
}
