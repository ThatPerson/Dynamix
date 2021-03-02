/**
 * @file datatypes.c
 */

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <time.h>    
//#define LOGGING

time_t start_time;

#ifdef LOGGING
#define LOG(a, args...) printf("LOG   %s(%s:%d %d) " a "\n",  __func__,__FILE__, __LINE__, (int) (time(0) - start_time), ##args)
#else
#define LOG(a, args...)
#endif

#define ERROR(a, args...) printf("ERROR %s(%s:%d %d) " a "\n",  __func__,__FILE__, __LINE__, (int) (time(0) - start_time), ##args)

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
#define MOD_AIMF 		23
#define MOD_AIMFT 		24

/*
 * For GAF/EGAF models. If VARIANT_A flag set in model struct, then the orientations of the orientation vectors
 * are allowed to vary.
 */
#define VARIANT_A		0
#define INVARIANT_A		1			


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
#define N_OR			14
#define OR_LIMIT		3.14

#define CNRATIO_ON		1
#define CNRATIO_OFF		0

#define N_RELAXATION	150						///< Approximate number of relaxation measurements; will dynamically allocate if overflows.
#define NTHREADS		40						///< Unused
#define THREAD_STACK	32768*2 				///< Bytes per thread. Raise if stack overflows, lower if insufficient stack for number of workers

#define GLOBAL			1
#define LOCAL			0

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

#define RDC_ON			1
#define RDC_OFF			0

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
	unsigned int params;
	struct Residue * residues;
	unsigned int n_residues;
	unsigned int rdc;
	unsigned int nthreads;
	unsigned int global;
	unsigned int cn_ratio;
	unsigned int or_variation; // VARIANT_A or INVARIANT_A.
	int error_mode;
	int proc_start;
	int proc_end;
	int myid, numprocs;
	double WS2NH, WS2CH, WS2CN, WS2CC;
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
 * @var Residue::S2XY
 *  Dipolar order parameter between X and Y
 * @var Residue::S2XYe
 *  Error in dipolar order parameter, twice standard deviation.
 * @var Residue::S2XYb
 *  Backup copy of order parameter for error calculation.
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
	double S2NHb, S2CCb, S2CHb, S2CNb; // order parameter backups for error calculation.
	double csisoN;
	double csisoC;
	double csaN[3];
	double csaC[3];
	struct Orient orients[14];
	//float orients[14][2]; // theta,phi
	struct Relaxation * relaxation;
	struct Relaxation * temp_relaxation;
	float cn;
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
	struct Model * model;
	unsigned int n_iter;
	char outputdir[255];
};


long double sq(long double x) {
	return x * x;
}

int sq_i(int x) {
	return x * x;
}

/**
 * Populates the second order spherical harmonics for orientation vector passed as pointer.
 * @param or
 *  Orientation vector
 * @return NULL
 */
void calculate_Y2(struct Orient * or) {
	/* or->Y2[0] is Y2 -2
	   or->Y2[1] is Y2 -1
	   or->Y2[2] is Y2 0
	   or->Y2[3] is Y2 1
	   or->Y2[4] is Y2 2 */

	/* as given https://mathworld.wolfram.com/SphericalHarmonic.html*/
	/* Y2[0] is Y2 -2*/

	or->Y2[4] = (1/4.) * (sqrt(15. / (2 * M_PI))) * (pow(sin(or->theta), 2.)) * cexp(2 * I * or->phi);
	or->Y2[3] = (-1/2.) * (sqrt(15. / (2 * M_PI))) * sin(or->theta) * cos(or->theta) * cexp(I * or->phi);
	or->Y2[2] = (1/4.) * (sqrt(5. / M_PI)) * (3 * pow(cos(or->theta), 2) - 1);
	or->Y2[1] = (1/2.) * (sqrt(15. / (2 * M_PI))) * sin(or->theta) * cos(or->theta) * cexp(-I * or->phi);
	or->Y2[0] = (1/4.) * (sqrt(15. / (2 * M_PI))) * (pow(sin(or->theta), 2.)) * cexp(-2 * I * or->phi);

	or->Y2c[0] = conj(or->Y2[0]);
	or->Y2c[1] = conj(or->Y2[1]);
	or->Y2c[2] = conj(or->Y2[2]);
	or->Y2c[3] = conj(or->Y2[3]);
	or->Y2c[4] = conj(or->Y2[4]);
	return;

}

/**
 * Initialises Wigner D matrix in global space. 
 */
void initialise_dwig(double angle, long double Dw[5][5]) {
	// beta taken as pi/2.
	/* It could be taken as an argument but hopefully
	 * doing it this way will mean the compiler will fill it out
	 * before execution, saving a little time. */
	// 0 -2
	// 1 -1
	// 2  0
	// 3  1
	// 4  2

	// d[m' + 2][m + 2]

	/* verified against https://link.springer.com/content/pdf/bbm%3A978-1-4684-0208-7%2F1.pdf in rotation_tests.c*/

	long double cosp = cosl(angle);
	long double sinp = sinl(angle);
	Dw[0][0] = powl(cosl(angle/2.), 4); // -2 -2 
	Dw[1][0] = (-1/2.) * (1 + cosp) * sinp; // -1 -2
	Dw[2][0] = sqrtl(3/8.) * powl(sinp, 2); // 0 -2
	Dw[3][0] = (1/2.) * sinp * (cosp - 1); // 1 -2
	Dw[4][0] = powl(sinl(angle/2.), 4); // 2 -2

	Dw[0][1] = (1/2.) * (1 + cosp) * sinp; // -2 -1
	Dw[1][1] = powl(cosp, 2.) - (1/2.) * (1 - cosp); // -1 -1
	Dw[2][1] = -sqrtl(3/8.) * sinl(angle * 2.); // 0 -1
	Dw[3][1] = (1/2.) * (1 + cosp) - powl(cosp, 2); // 1 -1
	Dw[4][1] = (1/2.) * (cosp - 1) * sinp; // 2 -1

	Dw[0][2] = sqrtl(3/8.) * powl(sinp, 2); // -2 0
	Dw[1][2] = sqrtl(3/2.) * sinp * cosp; // -1 0
	Dw[2][2] = (1/2.) * (3 * powl(cosp, 2) - 1); // 0 0
	Dw[3][2] = -sqrtl(3/2.) * sinp * cosp; // 1 0
	Dw[4][2] = sqrtl(3/8.) * powl(sinp, 2);  // 2 0

	Dw[0][3] = -(1/2.) * (cosp - 1) * sinp; // -2 1 
	Dw[1][3] = (1/2.) * (1 + cosp) - powl(cosp, 2); // -1 1
	Dw[2][3] = sqrtl(3/8.) * sinl(angle * 2.); // 0 1
	Dw[3][3] = powl(cosp, 2) - (1/2.)*(1-cosp); // 1 1
	Dw[4][3] = (-1/2.) * (1 + cosp) * sinp; // 2 1

	Dw[0][4] = powl(sinl(angle/2.), 4); // -2 2
	Dw[1][4] = (-1/2.) * (cosp - 1) * sinp; // -1 2
	Dw[2][4] = sqrtl(3/8.) * powl(sinp, 2); // 0 2
	Dw[3][4] = (1/2.) * sinp * (cosp + 1); // 1 2  m' = 1, m = 2
	Dw[4][4] = powl(cosl(angle / 2.), 4); // 2 2
	return;
}

/**
 * Function to free all allocated memory within Model.
 * @params m
 *  Pointer to Model to be freed.
 */
void free_all(struct Model *m) {
	unsigned int res, k;
	unsigned int params = m->params;

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

/**
 * Function to rotate spherical harmonics using Wigner D Matrices. I believe the
 * angles alpha, beta, gamma are aligned with the alpha, beta and gamma axis for deflections
 * (beta definitely is, I'm unsure about the others - but I'd assume so?)
 * @params or
 *  Pointer to Orient struct to modify.
 * @params alpha
 *  Alpha rotation
 * @params beta
 *  Beta rotation
 * @params gamma
 *  Gamma rotation
 * @returns Nothing
 *  Modifies orient inline.
 */
void rotate_Y2(struct Orient * or, double alpha, double beta, double gamma) {
	// the idea should be to run calculate_Y2 _first_, and then this function.
	// it will take the orient and overwrite it with the rotated vectors.

	long double Dw[5][5];

	initialise_dwig(beta, Dw); // we initialise Dw[][] with the reduced Wigner matrices.


	
	/* the Wigner D-matrix D^j_{m', m} is e^(-i m' alpha) d^j_{m', m} e^(-i m gamma) */

	double complex Y2b[5]; // we store Y2.
	int i;
	double complex alpha_m[5]; // precalculation of the exp(-I * mp * alpha) and exp(-I * m * gamma) parts.
	double complex gamma_m[5];
	for (i = 0; i < 5; i++) {
		Y2b[i] = or->Y2[i];
		or->Y2[i] = 0;
		or->Y2c[i] = 0;
		alpha_m[i] = cexp(-I * (i - 2) * alpha);
		gamma_m[i] = cexp(-I * (i - 2) * gamma);
	}

	int m, mp;
	double complex Dcomp = 0;

	for (mp = -2; mp <= 2; mp++) {
		for (m = -2; m <= 2; m++) {
			Dcomp = gamma_m[m+2];
			Dcomp *= Dw[mp+2][m+2];
			Dcomp *= alpha_m[mp+2]; 
			// D^j_{mp, m} = exp(-I * mp * alpha) * d^j_{m'm} * exp(-I * m * gamma)
			or->Y2[mp+2] += conj(Dcomp) * Y2b[m+2];
		}
		or->Y2c[mp+2] = conj(or->Y2[mp+2]); // and fill conjugate
	}

	return;
}
