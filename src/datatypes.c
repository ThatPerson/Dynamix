#include <stdlib.h>
#include <math.h>
#include <complex.h>
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


#define N_RELAXATION	50
#define NTHREADS		40	
#define THREAD_STACK	32768*2



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

#define MODE_REAL		0
#define MODE_IMAG		1
#define MODE_COMP		2

#define HALF_PI			M_PI / 2.

long double Dwig[5][5];



struct Model {
	int max_func_evals, max_iter, n_iter, n_error_iter;
	char outputdir[255];
	int model;
	struct Residue * residues;
	int n_residues;
	int nthreads;
	int error_mode;
};

struct Orient {
	float phi;
	float theta;
	long double complex Y2[5];
	long double complex Y2c[5];
};

struct Residue {
	float S2_dipolar;
	float csisoN;
	float csisoC;
	float csaN[3];
	float csaC[3];
	struct Orient orients[14];
	//float orients[14][2]; // theta,phi
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

void calculate_Y2(struct Orient * or) {
	or->Y2[0] = (1/4.) * (sqrtl(15. / (2 * M_PI))) * (powl(sinl(or->theta), 2.)) * cexpl(2 * I * or->phi);
	or->Y2[1] = (-1/2.) * (sqrtl(15. / (2 * M_PI))) * sinl(or->theta) * cosl(or->theta) * cexpl(I * or->phi);
	or->Y2[2] = (1/4.) * (sqrtl(5. / M_PI)) * (3 * powl(cosl(or->theta), 2) - 1);
	or->Y2[3] = (1/2.) * (sqrtl(15. / (2 * M_PI))) * sinl(or->theta) * cosl(or->theta) * cexpl(-I * or->phi);
	or->Y2[4] = (1/4.) * (sqrtl(15. / (2 * M_PI))) * (powl(sinl(or->theta), 2.)) * cexpl(-2 * I * or->phi);
	
	or->Y2c[0] = conjl(or->Y2[0]);
	or->Y2c[1] = conjl(or->Y2[1]);
	or->Y2c[2] = conjl(or->Y2[2]);
	or->Y2c[3] = conjl(or->Y2[3]);
	or->Y2c[4] = conjl(or->Y2[4]);
	
	return;
	
}


void initialise_dwig(void) {
	// beta taken as pi/2. 
	/* It could be taken as an argument but hopefully
	 * doing it this way will mean the compiler will fill it out
	 * before execution, saving a little time. */
	// 0 -2
	// 1 -1 
	// 2  0
	// 3  1
	// 4  2
	long double cosp = cosl(HALF_PI);
	long double sinp = sinl(HALF_PI);
	Dwig[0][0] = powl(cosl(HALF_PI/2.), 4);
	Dwig[1][0] = (-1/2.) * (1 + cosp) * sinp;
	Dwig[2][0] = sqrtl(3/8.) * powl(sinp, 2);
	Dwig[3][0] = (1/2.) * sinp * (cosp - 1);
	Dwig[4][0] = powl(sinl(HALF_PI/2.), 4);
	
	Dwig[0][1] = (1/2.) * (1 + cosp) * sinp;
	Dwig[1][1] = powl(cosp, 2.) - (1/2.) * (1 - cosp);
	Dwig[2][1] = -sqrtl(3/8.) * sinl(HALF_PI * 2.);
	Dwig[3][1] = (1/2.) * (1 + cosp) - powl(cosp, 2);
	Dwig[4][1] = (1/2.) * (cosp - 1) * sinp;
	
	Dwig[0][2] = sqrtl(3/8.) * powl(sinp, 2);
	Dwig[1][2] = sqrtl(3/2.) * sinp * cosp;
	Dwig[2][2] = (1/2.) * (3 * powl(cosp, 2) - 1);
	Dwig[3][2] = -sqrtl(3/2.) * sinp * cosp;
	Dwig[4][2] = sqrtl(3/8.) * powl(sinp, 2);
	
	Dwig[0][3] = -(1/2.) * (cosp - 1) * sinp;
	Dwig[1][3] = (1/2.) * (1 + cosp) - powl(cosp, 2);
	Dwig[2][3] = sqrtl(3/8.) * sinl(HALF_PI * 2.);
	Dwig[3][3] = powl(cosp, 2) - (1/2.)*(1-cosp);
	Dwig[4][3] = (-1/2.) * (1 + cosp) * sinp;
	
	Dwig[0][4] = powl(sinl(HALF_PI/2.), 4);
	Dwig[1][4] = (-1/2.) * (cosp - 1) * sinp;
	Dwig[2][4] = sqrtl(3/8.) * powl(sinp, 2);
	Dwig[3][4] = (1/2.) * sinp * (cosp + 1);
	Dwig[4][4] = powl(cosl(HALF_PI / 2.), 4);
	return;
}

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
		if (m->residues[res].relaxation != NULL)
			free(m->residues[res].relaxation);
		if (m->residues[res].parameters != NULL)
			free(m->residues[res].parameters);
		if (m->error_mode == 1) {
			if (m->residues[res].error_params != NULL) {
				for (k = 0; k < params; k++) {
					if (m->residues[res].error_params[k] != NULL)
						free(m->residues[res].error_params[k]);
				}
				free(m->residues[res].error_params);
			}
			if (m->residues[res].errors_mean != NULL)
				free(m->residues[res].errors_mean);
			if (m->residues[res].errors_std != NULL)
				free(m->residues[res].errors_std);
		}
	}
	free(m->residues);
	
	return;
}
