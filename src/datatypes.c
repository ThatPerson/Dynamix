/**
 * @file datatypes.c
 */

#include <stdlib.h>
#include <math.h>
#include <complex.h>
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
#define THREAD_STACK	(32768*2) 				///< Bytes per thread. Raise if stack overflows, lower if insufficient stack for number of workers

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


#define HALF_PI			(M_PI / 2.)

#define VERBOSE			1

#define T_DOWN			0.0000000010000
#define T_UP			1000000000
#define T_S				9

#define RDC_ON			1
#define RDC_OFF			0

Decimal Dwig[5][5]; 						///< 5x5 array containing Wigner components for pi/2


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
	unsigned int n_nm_iter, n_anneal_iter, n_error_iter;
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
	unsigned int proc_start;
	unsigned int proc_end;
	int myid, numprocs;
	Decimal WS2NH, WS2CH, WS2CN, WS2CC;

	double anneal_temp;
	double anneal_wobb;
	double anneal_therm;
	double anneal_restart;
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
	Decimal phi;
	Decimal theta;
	Complex Y2[5];
	Complex Y2c[5];
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
	Decimal S2NH, S2CC, S2CH, S2CN;
	Decimal S2NHe,S2CCe,S2CHe,S2CNe;
	Decimal S2NHb, S2CCb, S2CHb, S2CNb; // order parameter backups for error calculation.
	Decimal csisoN;
	Decimal csisoC;
	Decimal csaN[3];
	Decimal csaC[3];
	struct Orient orients[14];
	//Decimal orients[14][2]; // theta,phi
	struct Relaxation * relaxation;
	struct Relaxation * temp_relaxation;
	Decimal cn;
	Decimal ** error_params;

	unsigned int n_relaxation;
	unsigned int lim_relaxation;
	int ignore;
	Decimal min_val;
	Decimal * parameters;
	Decimal * errors_std;
	Decimal * errors_mean;
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
	Decimal field; // in MHz
	Decimal wr; // in Hz
	Decimal w1; // in Hz
	int type;
	Decimal R;
	Decimal Rerror;
	Decimal T; // in Kelvin
};

Decimal sq(Decimal x) {
	return x * x;
}

int sq_i(int x) {
	return x * x;
}

Decimal temp_tau(Decimal tau0, Decimal Ea, Decimal temp) {
    return tau0 * exp(Ea / (RYD * temp));
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

}

/**
 * Initialises Wigner D matrix in global space. 
 */
void initialise_dwig(Decimal angle, Decimal Dw[5][5]) {
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

	Decimal cosp = cos(angle);
	Decimal sinp = sin(angle);
	Dw[0][0] = pow(cos(angle/2.), 4); // -2 -2
	Dw[1][0] = (-1/2.) * (1 + cosp) * sinp; // -1 -2
	Dw[2][0] = sqrt(3/8.) * pow(sinp, 2); // 0 -2
	Dw[3][0] = (1/2.) * sinp * (cosp - 1); // 1 -2
	Dw[4][0] = pow(sin(angle/2.), 4); // 2 -2

	Dw[0][1] = (1/2.) * (1 + cosp) * sinp; // -2 -1
	Dw[1][1] = pow(cosp, 2.) - (1/2.) * (1 - cosp); // -1 -1
	Dw[2][1] = -sqrt(3/8.) * sin(angle * 2.); // 0 -1
	Dw[3][1] = (1/2.) * (1 + cosp) - pow(cosp, 2); // 1 -1
	Dw[4][1] = (1/2.) * (cosp - 1) * sinp; // 2 -1

	Dw[0][2] = sqrt(3/8.) * pow(sinp, 2); // -2 0
	Dw[1][2] = sqrt(3/2.) * sinp * cosp; // -1 0
	Dw[2][2] = (1/2.) * (3 * pow(cosp, 2) - 1); // 0 0
	Dw[3][2] = -sqrt(3/2.) * sinp * cosp; // 1 0
	Dw[4][2] = sqrt(3/8.) * pow(sinp, 2);  // 2 0

	Dw[0][3] = -(1/2.) * (cosp - 1) * sinp; // -2 1 
	Dw[1][3] = (1/2.) * (1 + cosp) - pow(cosp, 2); // -1 1
	Dw[2][3] = sqrt(3/8.) * sin(angle * 2.); // 0 1
	Dw[3][3] = pow(cosp, 2) - (1/2.)*(1-cosp); // 1 1
	Dw[4][3] = (-1/2.) * (1 + cosp) * sinp; // 2 1

	Dw[0][4] = pow(sin(angle/2.), 4); // -2 2
	Dw[1][4] = (-1/2.) * (cosp - 1) * sinp; // -1 2
	Dw[2][4] = sqrt(3/8.) * pow(sinp, 2); // 0 2
	Dw[3][4] = (1/2.) * sinp * (cosp + 1); // 1 2  m' = 1, m = 2
	Dw[4][4] = pow(cos(angle / 2.), 4); // 2 2
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
void rotate_Y2(struct Orient * or, Decimal alpha, Decimal beta, Decimal gamma) {
	// the idea should be to run calculate_Y2 _first_, and then this function.
	// it will take the orient and overwrite it with the rotated vectors.

	Decimal Dw[5][5];

	initialise_dwig(beta, Dw); // we initialise Dw[][] with the reduced Wigner matrices.


	
	/* the Wigner D-matrix D^j_{m', m} is e^(-i m' alpha) d^j_{m', m} e^(-i m gamma) */

	Complex Y2b[5]; // we store Y2.
	int i;
	Complex alpha_m[5]; // precalculation of the exp(-I * mp * alpha) and exp(-I * m * gamma) parts.
	Complex gamma_m[5];
	for (i = 0; i < 5; i++) {
		Y2b[i] = or->Y2[i];
		or->Y2[i] = 0;
		or->Y2c[i] = 0;
		alpha_m[i] = cexp(-I * (i - 2) * alpha);
		gamma_m[i] = cexp(-I * (i - 2) * gamma);
	}

	int m, mp;
	Complex Dcomp;

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

}

/**
 * Generates uniform random Decimal from 0 to 1 inclusive.
 * @return Decimal
 *  Long Decimal containing uniform random number.
 */
Decimal uniform_rand(void) {
    return ((Decimal) rand() + 1.) / ((Decimal) RAND_MAX + 1.);
}

void gen_params(const Decimal *minv, const Decimal *maxv, Decimal *pars, unsigned int n_pars) {
    unsigned int k;
    for (k = 0; k < n_pars; k++)
        pars[k] = minv[k] + (uniform_rand() * (maxv[k] - minv[k]));
}

void setup_paramlims(struct Model *m, Decimal S2NH, Decimal * minv, Decimal * maxv) {
    unsigned int k;
    for (k = 0; k < m->params; k++)
        minv[k] = 0;

    switch (m->model) {
        case MOD_SMF:
            maxv[0] = 10; maxv[1] = 1;
            break;
        case MOD_SMFT:
            maxv[0] = 1; maxv[1] = 1; maxv[2] = 60000;
            break;
        case MOD_EMF:
            minv[0] = 0.1; minv[1] = S2NH;
            maxv[0] = 10; maxv[1] = 1; maxv[2] = 0.1;
            break;
        case MOD_EMFT:
            minv[0] = pow(10, -7); minv[1] = S2NH;
            maxv[0] = pow(10, -4); maxv[1] = 1; maxv[2] = pow(10, -7); maxv[3] = 60000, maxv[4] = 60000;
            break;
        case MOD_DEMF:
            minv[0] = 0.1; minv[1] = S2NH; minv[3] = S2NH;
            maxv[0] = 10; maxv[1] = 1; maxv[2] = 0.1; maxv[3] = 1;
            break;
        case MOD_DEMFT:
            minv[0] = pow(10, -7); minv[1] = S2NH; minv[3] = S2NH;
            maxv[0] = pow(10, -4); maxv[1] = 1; maxv[2] = pow(10, -7); maxv[3] = 1; maxv[4] = 60000, maxv[5] = 60000;
            break;
        case MOD_GAF:
            minv[0] = 0.1;
            maxv[0] = 10; maxv[1] = 1;
            for (k = 2; k <= 7; k++) maxv[k] = 0.25;
            break;
        case MOD_GAFT:
            minv[0] = pow(10, -7);
            maxv[0] = pow(10, -4); maxv[1] = pow(10, -7);
            for (k = 2; k <= 7; k++) maxv[k] = 0.25;
            maxv[8] = 60000; maxv[9] = 60000;
            break;
        case MOD_AIMF:
            minv[0] = 0.1;
            maxv[0] = 10; maxv[1] = 1;
            for (k = 2; k <= 7; k++) { maxv[k] = 1; minv[k] = S2NH; }
            break;
        case MOD_AIMFT:
            minv[0] = pow(10, -7);
            maxv[0] = pow(10, -4); maxv[1] = pow(10, -7);
            for (k = 2; k <= 7; k++) { maxv[k] = 1; minv[k] = S2NH; }
            maxv[8] = 60000; maxv[9] = 60000;
            break;
        case MOD_EGAF:
            minv[0] = 0.1;
            maxv[0] = 10; maxv[1] = 1;
            for (k = 2; k <= 4; k++) maxv[k] = 0.25;
            minv[5] = S2NH; maxv[5] = 1;
            break;
        case MOD_EGAFT:
            minv[0] = pow(10, -7);
            maxv[0] = pow(10, -4); maxv[1] = pow(10, -7);
            for (k = 2; k <= 4; k++) maxv[k] = 0.25;
            minv[5] = S2NH; maxv[5] = 1;
            maxv[6] = 60000; maxv[7] = 60000;
            break;
        default: break;
    }
    if (m->or_variation == VARIANT_A && m->rdc == RDC_ON) {
        /* eg for MOD_GAF, params = 8 + 3. opts[7] is full, so we want to put
         * alpha in opts[params-3] and beta in opts[params-2] and gamma in opts[params-1];
         */
        minv[m->params - 6] = 0; // papbC
        minv[m->params - 5] = 0; // papbN
        minv[m->params - 4] = (rand() % 20000) / 1.; // kex
        minv[m->params - 3] = 0; // alpha
        minv[m->params - 2] = 0; // beta
        minv[m->params - 1] = 0; // gamma
        maxv[m->params - 6] = 0; // papbC
        maxv[m->params - 5] = 0; // papbN
        maxv[m->params - 4] = (rand() % 20000) / 1.; // kex
        maxv[m->params - 3] = 0; // alpha
        maxv[m->params - 2] = 0; // beta
        maxv[m->params - 1] = 0; // gamma
    } else if (m->or_variation == VARIANT_A) {
        minv[m->params - 3] = 0; // alpha
        minv[m->params - 2] = 0; // beta
        minv[m->params - 1] = 0; // gamma
        maxv[m->params - 3] = 0; // alpha
        maxv[m->params - 2] = 0; // beta
        maxv[m->params - 1] = 0; // gamma
    } else if (m->rdc == RDC_ON) {
        minv[m->params - 3] = 0; // papbC
        minv[m->params - 2] = 0; // papbN
        minv[m->params - 1] = (rand() % 20000) / 1.; // kex#
        maxv[m->params - 3] = 0; // papbC
        maxv[m->params - 2] = 0; // papbN
        maxv[m->params - 1] = (rand() % 20000) / 1.; // kex
    }
/**
         * SMF parameters are \n
         *   [0] tau\n
         *   [1] S2\n
         * SMFT parameters;\n
         *   [0] tau\n
         *   [1] S2\n
         *   [2] Ea\n
         * EMF parameters\n
         *   [0] tau slow\n
         *   [1] S2 slow\n
         *   [2] tau fast\n
         *   NOTE: The fast order parameter is calculated as S2NH/S2s\n
         * EMFT parameters\n
         *   [0] tau slow\n
         *   [1] S2 slow\n
         *   [2] tau fast\n
         *   [3] activation energy for slow motion\n
         *   [4] activation energy for fast motion\n
         * DEMF parameters\n
         *   [0] tau slow\n
         *   [1] S2 slow\n
         *   [2] tau fast\n
         *   [3] S2 fast\n
         * DEMFT parameters\n
         *   [0] tau slow\n
         *   [1] S2 slow\n
         *   [2] tau fast\n
         *   [3] S2 fast\n
         *   [4] activation energy for slow motion\n
         *   [5] activation energy for fast motion\n
         * GAF parameters\n
         *   [0] tau slow\n
         *   [1] tau fast\n
         *   [2-4] alpha, beta, gamma deflections for slow motions\n
         *   [5-7] alpha, beta, gamma deflections for fast motions\n
         * GAFT parameters\n
         *   [0] tau slow\n
         *   [1] tau fast\n
         *   [2-4] alpha, beta, gamma deflections for slow motions\n
         *   [5-7] alpha, beta, gamma deflections for fast motions\n
         *   [8] activation energy for slow motion\n
         *   [9] activation energy for fast motion\n
         * EGAF parameters\n
         *   [0] tau slow\n
         *   [1] tau fast\n
         *   [2-4] alpha, beta, gamma deflections for slow motions\n
         *   [5] fast motion order parameter\n
         * EGAFT parameters\n
         *   [0] tau slow\n
         *   [1] tau fast\n
         *   [2-4] alpha, beta, gamma deflections for slow motions\n
         *   [5] order parameter for fast motions\n
         *   [6] activation energy for slow motion\n
         *   [7] activation energy for fast motion\n
         */
}

