#ifndef DATATYPES_H
#define DATATYPES_H

#include <time.h>
#include <math.h>

typedef double Decimal;
typedef double _Complex Complex;

#define MPI_DECIMAL MPI_DOUBLE

#ifdef LOGGING
#define LOG(a, args...) printf("LOG   %s(%s:%d %d) " a "\n",  __func__,__FILE__, __LINE__, (int) (time(0) - start_time), ##args)
#else
#define LOG(a, args...) /* nothing */
#endif

#define ERROR(a, args...) printf("ERROR %s(%s:%d %d) " a "\n",  __func__,__FILE__, __LINE__, (int) (time(0) - start_time), ##args)

#define MOD_UNDEFINED    0
#define MOD_SMF        1                        ///< Simple Model Free
#define MOD_EMF        2                        ///< Extended Model Free using S2 dipolar
#define MOD_EMFT        3                        ///< Extended Model Free with Temperature Dependence using S2 dipolar
#define MOD_SMFT        4                        ///< Simple Model Free with Temperature Dependence
#define MOD_DEMF        5                        ///< Extended Model Free
#define MOD_DEMFT        6                        ///< Extended Model Free  with Temperature Dependence
#define MOD_GAF        7                        ///< Gaussian Axial Fluctuations model
#define MOD_GAFT        8                        ///< Gaussian Axial Fluctuations model with Temperature Dependence
#define MOD_EGAF        9                        ///< 3D GAF for slow motions, EMF for fast motions
#define MOD_EGAFT        10                        ///< 3D GAF for slow motions, EMF for fast motions, temperature dependence.
#define MOD_AIMF        23
#define MOD_AIMFT        24
#define MOD_BGF         25
#define MOD_BGFT        26
#define MOD_EAIMF       27
#define MOD_EAIMFT      28
#define MOD_BGAF        29
#define MOD_BGAFT       30
#define MOD_BAIMF       31
#define MOD_BAIMFT      32
#define MOD_RDEMFT      33
/*
 * For GAF/EGAF models. If VARIANT_A flag set in model struct, then the orientations of the orientation vectors
 * are allowed to vary.
 */
#define VARIANT_A        0
#define INVARIANT_A        1

#define ENABLED 1
#define DISABLED 0

#define DATA_S2NH        0
#define DATA_S2CH        1
#define DATA_S2CC        2
#define DATA_S2CN        3
#define DATA_CSISON        4
#define DATA_CSISOC        5

#define R_15NR1        0
#define R_15NR1p        1
#define R_13CR1        2
#define R_13CR1p        3

#define OR_NH            0
#define OR_NC            1
#define OR_NCA            2
#define OR_NCSAxx        3
#define OR_NCSAyy        4
#define OR_NCSAzz        5
#define OR_CCAp        6
#define OR_CCAc        7
#define OR_CN            8
#define OR_CNH            9
#define OR_CH            10
#define OR_CCSAxx        11
#define OR_CCSAyy        12
#define OR_CCSAzz        13
#define N_OR            14
#define OR_LIMIT        3.14

#define CNRATIO_ON        1
#define CNRATIO_OFF        0

#define N_RELAXATION    150                        ///< Approximate number of relaxation measurements; will dynamically allocate if overflows.

#define GLOBAL            1
#define LOCAL            0

#define MIN_VAL            10000000                ///< Maximum value of chisq function for which it will be considered fit

#define RYD            8.3144621                ///< Rydberg Constant in JK-1mol-1

/* Dipolar couplings are calculated using
 * D = (10^-7 * gammaX * gammaY * (h/(2pi))) / (r^3)
 * where gammaX, gammaY are in units of rad / s T
 *   and r in units of Angstroms.
 * See utils/calc_D.py
 */
#define D_NH    72038.411068
#define D_NC    8175.249117
#define D_NCA    6180.126304
#define D_CCAp    13467.660521
#define D_CCAc    3093.881320
#define D_CN    8175.249117
#define D_CH    22362.477236
#define D_NHr    13108.322725
#define D_CHr    31491.710460

#define MODE_15N        0
#define MODE_13C        1

#define MODE_REAL        0
#define MODE_IMAG        1


#define HALF_PI            (M_PI / 2.)

#define VERBOSE            1

#define T_DOWN            0.0000000010000
#define T_UP            1000000000
#define T_S                9

#define DEUTERATED 0
#define PROTONATED 1

extern time_t start_time;
extern Decimal Dwig[5][5];
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
    struct Residue *residues;
    unsigned int n_residues;
    unsigned int nthreads;
    unsigned int global;
    unsigned int cn_ratio;
    unsigned int or_variation; // VARIANT_A or INVARIANT_A.
    int error_mode;
    unsigned int proc_start;
    unsigned int proc_end;
    int myid, numprocs;
    int WS2NH, WS2CH, WS2CN, WS2CC;
    unsigned int ultrafast;
    unsigned int microsecond;
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
    Decimal rot_phi, rot_theta;
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
    Decimal S2NHe, S2CCe, S2CHe, S2CNe;
    Decimal S2NHb, S2CCb, S2CHb, S2CNb; // order parameter backups for error calculation.
    Decimal csisoN;
    Decimal csisoC;
    Decimal csaN[3];
    Decimal csaC[3];
    struct Orient orients[14];
    //Decimal orients[14][2]; // theta,phi
    struct Relaxation *relaxation;
    struct Relaxation *temp_relaxation;
    Decimal cn;
    Decimal **error_params;

    unsigned int n_relaxation;
    unsigned int lim_relaxation;
    int ignore;
    Decimal min_val;
    Decimal *parameters;
    Decimal *errors_std;
    Decimal *errors_mean;
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
    int hydrogen;
};

Decimal temp_tau(Decimal tau0, Decimal Ea, Decimal temp);

void calculate_Y2(struct Orient *or);

void initialise_dwig(Decimal angle, Decimal Dw[5][5]);

void free_all(struct Model *m);

Decimal sq(Decimal x);

void rotate_Y2(struct Orient *or, Decimal alpha, Decimal beta, Decimal gamma);

Decimal uniform_rand(void);

void gen_params(const Decimal *minv, const Decimal *maxv, Decimal *pars, unsigned int n_pars);

void setup_paramlims(struct Model *m, Decimal S2NH, Decimal *minv, Decimal *maxv);

int determine_residues(unsigned int n_res, int myid, int numprocs, unsigned int *start, unsigned int *end);

#endif
