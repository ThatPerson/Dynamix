/**
 * @file datatypes.c
 */

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define MOD_UNDEFINED	0
#define MOD_SMF 		1						///< Simple Model Free
#define MOD_EMF 		2						///< Extended Model Free using S2 dipolar
#define MOD_EMFT		3						///< Extended Model Free with Temperature Dependence using S2 dipolar
#define MOD_SMFT 		4						///< Simple Model Free with Temperature Dependence
#define MOD_DEMF		5 						///< Extended Model Free 
#define MOD_DEMFT		6						///< Extended Model Free  with Temperature Dependence
#define MOD_GAF 		7						///< Gaussian Axial Fluctuations model
#define MOD_GAFT 		8						///< Gaussian Axial Fluctuations model with Temperature Dependence
#define MOD_EGAF		9						///< 3D GAF for slow motions, EMF for fast motions
#define MOD_EGAFT		10						///< 3D GAF for slow motions, EMF for fast motions, temperature dependence.

#define DATA_S2NH		0
#define DATA_S2CH		1
#define DATA_S2CC		2
#define DATA_S2CN		3
#define DATA_CSISON		4
#define DATA_CSISOC		5

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
#define OR_CCAp		6
#define OR_CCAc		7
#define OR_CN			8
#define OR_CNH			9
#define OR_CH			10
#define OR_CCSAxx		11
#define OR_CCSAyy		12
#define OR_CCSAzz		13


#define N_RELAXATION	150						///< Approximate number of relaxation measurements; will dynamically allocate if overflows.
#define NTHREADS		40						///< Unused
#define THREAD_STACK	32768*2 				///< Bytes per thread. Raise if stack overflows, lower if insufficient stack for number of workers



#define MIN_VAL			10000000				///< Maximum value of chisq function for which it will be considered fit

#define RYD 			8.3144621 				///< Rydberg Constant in JK-1mol-1

/* These constants are extracted from NHCO_3GAFSquaredEa, using the bond lengths defined in there */
#define D_NH			72038.41107 			///< Dipolar couplings taken from NHCO_3GAFSquaredEa for amide N-H
#define D_CC			13467.66052 			///< Dipolar couplings taken from NHCO_3GAFSquaredEa for alphatic C-Calpha
#define D_CH			22362.47724 			///< Dipolar couplings taken from NHCO_3GAFSquaredEa for C-H
#define D_CHr			31491.71046 			///< Dipolar couplings taken from NHCO_3GAFSquaredEa for C-Hrest
#define D_CN			8175.2 					///< Dipolar couplings taken from NHCO_3GAFSquaredEa for C-N
#define D_CaN			6180.1 					///< Dipolar couplings taken from NHCO_3GAFSquaredEa for amide N - aliphatic C
#define D_HNr			13108.32273 			///< Dipolar couplings taken from NHCO_3GAFSquaredEa for amide N-Hrest

#define MODE_15N		0
#define MODE_13C		1

#define MODE_REAL		0
#define MODE_IMAG		1
#define MODE_COMP		2

#define DIPOLAR			1
#define NODIP			0

#define HALF_PI			M_PI / 2.

#define VERBOSE			1

#define T_DOWN			0.0000000010000
#define T_UP			1000000000
#define T_S				9


long double Dwig[5][5]; 						///< 5x5 array containing Wigner components for pi/2


/** @struct Model
   *  Struct containing model parameters
   *
   *  @var Model::max_func_evals
   *    Unused
   *  @var Model::max_iter
   *    Unused
   *  @var Model::n_iter
   *    Number of iterations in main channel
   *  @var Model::n_error_iter
   *    Number of iterations for error calculation
   *  @var Model::outputdir
   *    Output directory
   *  @var Model::model
   *    Defined as one of MOD_*
   *  @var Model::residues
   *    Pointer to an array of n_residues residues
   *  @var Model::n_residues
   *    Length of residues; number of residues in protein
   *  @var Model::nthreads
   *    Number of threads to spawn for calculations
   *  @var Model::error_mode
   *    1 indicates that errors should be calculated, 0 indicates not.
   */
struct Model {
	unsigned int max_func_evals, max_iter, n_iter, n_error_iter;
	char outputdir[255];
	unsigned int model;
	struct Residue * residues;
	unsigned int n_residues;
	unsigned int nthreads;
	int error_mode;
};

/** @struct Orient
 * Struct defining atomic relative orientations
 * @var Orient::phi
 *  Phi angle (radians)
 * @var Orient::theta
 *  Theta angle (radians)
 * @var Orient::Y2
 *  Array containing the second order spherical harmonics, Y2m, for m = -2 (0) to 2 (5).
 * @var Orient::Y2c
 *  Conjugate to Y2.
 */
struct Orient {
	float phi;
	float theta;
	double complex Y2[5];
	double complex Y2c[5];
};

/** @struct Residue
 * Struct defining a residue
 * @var Residue::S2_dipolar
 *  Dipolar order parameter; used in EMF analysis
 * @var Residue::csisoN
 *  Isotropic chemical shift component for amide nitrogen.
 * @var Residue::csisoC
 *  Isotropic chemical shift component for carbonyl carbon.
 * @var Residue::csaN
 *  Anisotropic chemical shift components for nitrogen
 * @var Residue::csaC
 *  Anisotropic chemical shift components for carbon
 * @var Residue::orients
 *  Orientations for intraresidue atom vectors, defined 0-13 by OR_*
 * @var Residue::relaxation
 *  Dynamically allocated pointer to array of structs containing the relaxation data
 * @var Residue::temp_relaxation
 *  Pointer used to keep track of relaxation data during error calculation
 * @var Residue::error_params
 *  Pointer to array of pointers, each of which points to set of optimized parameters for each error iteration
 * @var Residue::n_relaxation
 *  Number of relaxation measurements for this residue
 * @var Residue::lim_relaxation
 *  Length of Residue::relaxation array
 * @var Residue::ignore
 *  Flag to ignore residue from calculations
 * @var Residue::min_val
 *  Optimized minimum chisq value
 * @var Residue::parameters
 *  Optimized parameters relating to minimum chisq value
 * @var Residue::errors_std
 *  Pointer to array containing standard deviations for each error parameter
 * @var Residue::errors_mean
 *  Pointer to array containing means for each error parameter
 * @var Residue::error_calcs
 *  Actual number of points used in error calculation
 */
struct Residue {
	double S2NH, S2CC, S2CH, S2CN;
	double S2NHe,S2CCe,S2CHe,S2CNe;
	double csisoN;
	double csisoC;
	double csaN[3];
	double csaC[3];
	struct Orient orients[14];
	//float orients[14][2]; // theta,phi
	struct Relaxation * relaxation;
	struct Relaxation * temp_relaxation;

	long double ** error_params;

	unsigned int n_relaxation;
	unsigned int lim_relaxation;
	int ignore;
	double min_val;
	long double * parameters;
	long double * errors_std;
	long double * errors_mean;
	int error_calcs;
};

/**
 * @struct Relaxation
 * Struct defining a relaxation measurement
 * @var Relaxation::field
 *  Proton nutation frequency for field in which relaxation measurement was made in MHz
 * @var Relaxation::wr
 *  Spinning frequency (MAS) for relaxation measurement (in Hz)
 * @var Relaxation::w1
 *  Nutation frequency for spin lock relaxation measurement (in Hz)
 * @var Relaxation::type
 *  One of R_15NR1, R_15NR1p, R_13CR1, R_13CR1p depending on which parameter this is a measurement of
 * @var Relaxation::R
 *  Relaxation measurement in s^-1
 * @var Relaxation::Rerror
 *  Monte-Carlo relaxaxtion measurement error (2 standard deviations) in s^-1
 * @var Relaxation::T
 *  Temperature of measurement (Kelvin)
 */
struct Relaxation {
	double field; // in MHz
	double wr; // in Hz
	double w1; // in Hz
	int type;
	double R;
	double Rerror;
	double T; // in Kelvin
};

/**
 * @struct rrargs
 * Struct containing information relevant to a thread
 * @var rrargs::i
 *  Count variable
 * @var rrargs::resid
 *  Pointer to residue being considered
 * @var rrargs::model
 *  Model type, one of MOD_*
 * @var rrargs::n_iter
 *  Number of iterations
 * @var rrargs::outputdir
 *  Output directory
 */
struct rrargs {
	unsigned int i;
	struct Residue * resid;
	unsigned int model;
	unsigned int n_iter;
	char outputdir[255];
};

/**
 * Populates the second order spherical harmonics for orientation vector passed as pointer.
 * @param or
 *  Orientation vector
 * @return NULL
 */
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

/**
 * Initialises Wigner D matrix in global space. 
 */
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

/**
 * Function to free all allocated memory within Model.
 * @params m
 *  Pointer to Model to be freed.
 */
void free_all(struct Model *m) {
	unsigned int res, k;
	unsigned int params;
	switch (m->model) {
		case MOD_SMF: params = 2; break;
		case MOD_EMF: params = 3; break;
		case MOD_EMFT:params = 5; break;
		case MOD_SMFT:params = 3; break;
		case MOD_DEMF: params = 4; break;
		case MOD_DEMFT: params= 6; break;
		case MOD_GAF: params = 8; break;
		case MOD_GAFT: params = 10; break;
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
