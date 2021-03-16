//
// Created by ben on 16/03/2021.
//

#include "model.h"
#include "../datatypes.h"

/** Calculates spectral density function for given frequency according to extended model free analysis.
 *  Implements J(w) = (1 - S2f) tauf / (1 + (w tauf)^2) + S2f (1 - S2s) taus / (1 + (w taus)^2),
 *  as in eq 2, Clore 1990 (assuming S2 = S2s*S2f and that there is no overall tumbling).
 */
Decimal J0(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f) {
    return ( \
			((((1 - (Decimal) S2f)) * (Decimal) tauf)) \
			/ (1 + ((Decimal) omega * (Decimal) tauf * (Decimal) omega * (Decimal) tauf)) \
		)\
		+\
		(\
			(((Decimal) S2f) * (1 - (Decimal) S2s) * (Decimal) taus)\
			/ (1 + ((Decimal) omega * (Decimal) taus * (Decimal) omega * (Decimal) taus))\
		);
}

/** Cross correlated spectral density. Taken from Manuscript,
 * Juv(w) = (1-S2f) tf / (1 + (wtf)^2) + (1/P2) S2f (P2 - S2s) (ts / 1+(wts)^2)
 * eq 2
 */
Decimal J0_CC(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f) {
    return ( \
			((((1 - (Decimal) S2f)) * (Decimal) tauf)) \
			/ (1 + ((Decimal) omega * (Decimal) tauf * (Decimal) omega * (Decimal) tauf)) \
		)\
		-\
		2. * (\
			(((Decimal) S2f) * (-0.5 - (Decimal) S2s) * (Decimal) taus)\
			/ (1 + ((Decimal) omega * (Decimal) taus * (Decimal) omega * (Decimal) taus))\
		);
}


Decimal Dipolar_R1(Decimal omega_obs, Decimal omega_neigh, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal D) {
    // for SMF and SMFT models, set S2s = 0 and taus = 0.
    Decimal q = (0.1) * sq(D) * ( \
		(J0(omega_neigh - omega_obs, taus, S2s, tauf, S2f)) \
		 + 3 * (J0(omega_obs, taus, S2s, tauf, S2f)) \
		  + 6 * (J0(omega_neigh + omega_obs, taus, S2s, tauf, S2f)) \
	  );
    return (Decimal) q;
}


Decimal Dipolar_R2(Decimal omega_obs, Decimal omega_neigh, Decimal w1, Decimal wr, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal D) {
    // for SMF and SMFT models, set S2s = 0 and taus = 0.
    return (Decimal) (\
  		(1/20.) * sq(D) * (\
  			(2/3.) * J0(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f) +\
  			(2/3.) * J0(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J0(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J0(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f) +\
  			(3.)   * J0(omega_obs, taus, S2s, tauf, S2f) +\
  			(1.)   * J0(omega_neigh - omega_obs, taus, S2s, tauf, S2f) +\
  			(6.)   * J0(omega_neigh, taus, S2s, tauf, S2f) +\
  			(6.)   * J0(omega_neigh + omega_obs, taus, S2s, tauf, S2f)\
  		)\
  	);
}

/**
 * Calculates the R2 contribution of CSA of nuclei with frequency omega.
 * D2 is the squared D22x/D22y/D22xy variable (dependent on CSA)
 * Kurbanov2011, eqs 18-20.
 */
Decimal CSA_R2(Decimal omega, \
						 Decimal w1, \
						 Decimal wr, \
						 Decimal taus, \
						 Decimal S2s, \
						 Decimal tauf, \
						 Decimal S2f, \
						 Decimal D2, \
						 Decimal  (*J_SD)(\
						 	Decimal,\
						 	Decimal, \
						 	Decimal, \
						 	Decimal, \
						 	Decimal)\
						) {
    return (Decimal) ( \
  		(1/45.) * (Decimal) D2 * ( \
  			(2/3.) * J_SD(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f) +\
  			(2/3.) * J_SD(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J_SD(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f) +\
  			(4/3.) * J_SD(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f) + \
  			(3.)   * J_SD(omega, taus, S2s, tauf, S2f)\
  		)\
  	);
}

struct BCParameters {
    Decimal S2NHs, S2NHf;
    Decimal S2NCSAxs, S2NCSAxf;
    Decimal S2NCSAys, S2NCSAyf;
    Decimal S2NCSAxys, S2NCSAxyf;
    Decimal S2CSAxs, S2CSAxf;
    Decimal S2CSAys, S2CSAyf;
    Decimal S2CSAxys, S2CSAxyf;
    Decimal S2CNs, S2CNf;
    Decimal S2CaNs, S2CaNf;
    Decimal taus, tauf;
    Decimal S2CHs, S2CHf;
    Decimal S2CCs, S2CCf;
};

Decimal Calc_15NR1(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars) {
    /* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
    Decimal field = relax->field * 1000000; // conversion to Hz

    /* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
     * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
     * software division) I've predefined it. Bond length taken as 1.02 A */

    Decimal omega_1H = 2 * M_PI * field;
    Decimal omega_13C, omega_15N;
    Decimal d2x, d2y, d2xy;

    omega_13C = 2 * M_PI * field / 3.976489314034722;
    omega_15N = 2 * M_PI * field / 9.869683408806043;

    Decimal *csa;
    csa = res->csaN;

    omega_1H *= T_DOWN;
    omega_13C *= T_DOWN;
    omega_15N *= T_DOWN;

    /* N CSA relaxation contribution */
    Decimal R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1NH, R1NHr, R1CN, R1CaN;
    Decimal J1;

    d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N * T_UP);
    d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N * T_UP);
    d2xy= (Decimal) sq(0.000001 * omega_15N * T_UP) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

    J1 = J0(omega_15N, pars->taus, pars->S2NCSAxs, pars->tauf, pars->S2NCSAxf);
    R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
    J1 = J0(omega_15N, pars->taus, pars->S2NCSAys, pars->tauf, pars->S2NCSAyf);
    R1CSAy = (1/15.) * d2y * J1;
    J1 = J0_CC(omega_15N, pars->taus, pars->S2NCSAxys, pars->tauf, pars->S2NCSAxyf);
    R1CSAxy = (1/15.) * d2xy * J1;
    /** Eq A30, Bremi1997 */
    R1CSA = R1CSAx + R1CSAy + 2.*R1CSAxy;

    /* N Dipolar Interactions Contributions */
    R1NH = Dipolar_R1(omega_15N, omega_1H, pars->taus, pars->S2NHs, pars->tauf, pars->S2NHf, D_NH);
    R1NHr = Dipolar_R1(omega_15N, omega_1H, pars->taus, pars->S2CNs, pars->tauf, pars->S2CNf, D_HNr);
    R1CN = Dipolar_R1(omega_15N, omega_13C, pars->taus, pars->S2CNs, pars->tauf, pars->S2CNf, D_CN);
    R1CaN = Dipolar_R1(omega_15N, omega_13C, pars->taus, pars->S2CaNs, pars->tauf, pars->S2CaNf, D_CaN);

    return (Decimal) (R1CSA + R1NH + R1NHr + R1CN + R1CaN) * T_DOWN;
}

Decimal Calc_15NR2(struct Residue *res, struct Relaxation* relax, struct BCParameters *pars) {
    /* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
    Decimal field = relax->field * 1000000; // conversion to Hz
    Decimal w1 = relax->w1, wr = relax->wr;
    /* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
     * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
     * software division) I've predefined it. Bond length taken as 1.02 A */

    Decimal omega_1H = 2 * M_PI * field;
    Decimal omega_13C, omega_15N;
    Decimal d2x, d2y, d2xy;

    omega_13C = 2 * M_PI * field / 3.976489314034722;
    omega_15N = 2 * M_PI * field / 9.869683408806043;

    Decimal *csa;
    csa = res->csaN;

    d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
    d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
    d2xy= (Decimal) sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

    /* N CSA relaxation contribution */
    Decimal R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2NH, R2NHr, R2CN, R2CaN;

    w1 *= T_DOWN;
    wr *= T_DOWN;
    omega_1H *= T_DOWN;
    omega_13C *= T_DOWN;
    omega_15N *= T_DOWN;

    R2CSAx = CSA_R2(omega_15N, w1, wr, pars->taus, pars->S2NCSAxs, pars->tauf, pars->S2NCSAxf, d2x, J0);
    R2CSAy = CSA_R2(omega_15N, w1, wr, pars->taus, pars->S2NCSAys, pars->tauf, pars->S2NCSAyf, d2y, J0);
    R2CSAxy = CSA_R2(omega_15N, w1, wr, pars->taus, pars->S2NCSAxys, pars->tauf, pars->S2NCSAxyf, d2xy, J0_CC);
    /**
     * Equation for spectral density for xy taken from Manuscript
     */
    R2CSA = R2CSAx + R2CSAy + 2.* R2CSAxy;

    /* N Dipolar Interactions Contributions */
    R2NH = Dipolar_R2(omega_15N, omega_1H, w1, wr, pars->taus, pars->S2NHs, pars->tauf, pars->S2NHf, D_NH);
    R2NHr = Dipolar_R2(omega_15N, omega_1H, w1, wr, pars->taus, pars->S2CNs, pars->tauf, pars->S2CNf, D_HNr);
    R2CN = Dipolar_R2(omega_15N, omega_13C, w1, wr, pars->taus, pars->S2CNs, pars->tauf, pars->S2CNf, D_CN);
    R2CaN = Dipolar_R2(omega_15N, omega_13C, w1, wr, pars->taus, pars->S2CaNs, pars->tauf, pars->S2CaNf, D_CaN);
    return (R2CSA + R2NH + R2NHr + R2CN + R2CaN)*T_DOWN;
}

Decimal Calc_13CR1(struct Residue *res, struct Relaxation* relax, struct BCParameters *pars) {
    /* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
    Decimal field = relax->field * 1000000; // conversion to Hz

    /* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
     * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
     * software division) I've predefined it. Bond length taken as 1.02 A */

    Decimal omega_1H = 2 * M_PI * field;
    Decimal omega_13C, omega_15N, wCOCa;
    Decimal d2x, d2y, d2xy;

    omega_13C = 2 * M_PI * field / 3.976489314034722;
    omega_15N = 2 * M_PI * field / 9.869683408806043;
    wCOCa = 120 * omega_13C * 0.000001;

    Decimal *csa;
    csa = res->csaC;

    d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
    d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
    d2xy= (Decimal) sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

    /* N CSA relaxation contribution */
    Decimal R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1CH, R1CHr, R1CN, R1CC;
    Decimal J1;


    omega_1H *= T_DOWN;
    omega_13C *= T_DOWN;
    omega_15N *= T_DOWN;

    J1 = J0(omega_13C, pars->taus, pars->S2CSAxs, pars->tauf, pars->S2CSAxf);
    R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
    J1 = J0(omega_13C, pars->taus, pars->S2CSAys, pars->tauf, pars->S2CSAyf);
    R1CSAy = (1/15.) * d2y * J1;
    J1 = J0_CC(omega_13C, pars->taus, pars->S2CSAxys, pars->tauf, pars->S2CSAxyf);
    R1CSAxy = (1/15.) * d2xy * J1;
    R1CSA = R1CSAx + R1CSAy + 2.*R1CSAxy;



    /* N Dipolar Interactions Contributions */
    R1CH = Dipolar_R1(omega_13C, omega_1H, pars->taus, pars->S2CHs, pars->tauf, pars->S2CHf, D_CH);
    R1CHr = Dipolar_R1(omega_13C, omega_1H, pars->taus, pars->S2CNs, pars->tauf, pars->S2CNf, D_CHr);
    R1CN = Dipolar_R1(omega_13C, omega_15N, pars->taus, pars->S2CNs, pars->tauf, pars->S2CNf, D_CN);
    R1CC = Dipolar_R1(omega_13C, omega_13C - wCOCa, pars->taus, pars->S2CCs, pars->tauf, pars->S2CCf, D_CC);

    return (Decimal) (R1CSA + R1CH + R1CHr + R1CN + R1CC)*T_DOWN;
}

Decimal Calc_13CR2(struct Residue *res, struct Relaxation* relax, struct BCParameters *pars) {
    /* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
    Decimal field = relax->field * 1000000; // conversion to Hz
    Decimal w1 = relax->w1, wr = relax->wr;
    /* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
     * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
     * software division) I've predefined it. Bond length taken as 1.02 A */

    Decimal omega_1H = 2 * M_PI * field;
    Decimal omega_13C, omega_15N, wCOCa;
    Decimal d2x, d2y, d2xy;

    omega_13C = 2 * M_PI * field / 3.976489314034722;
    omega_15N = 2 * M_PI * field / 9.869683408806043;
    wCOCa = 120 * omega_13C * 0.000001;

    Decimal *csa;
    csa = res->csaC;

    d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C);
    d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C);
    d2xy= (Decimal) sq(0.000001 * omega_13C) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

    /* CSA relaxation contribution */
    Decimal R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2CH, R2CHr, R2CN, R2CC;

    w1 *= T_DOWN;
    wr *= T_DOWN;
    omega_1H *= T_DOWN;
    omega_13C *= T_DOWN;
    omega_15N *= T_DOWN;

    R2CSAx = CSA_R2(omega_13C, w1, wr, pars->taus, pars->S2CSAxs, pars->tauf, pars->S2CSAxf, d2x, J0);
    R2CSAy = CSA_R2(omega_13C, w1, wr, pars->taus, pars->S2CSAys, pars->tauf, pars->S2CSAyf, d2y, J0);
    R2CSAxy = CSA_R2(omega_13C, w1, wr, pars->taus, pars->S2CSAxys, pars->tauf, pars->S2CSAxyf, d2xy, J0_CC);
    R2CSA = R2CSAx + R2CSAy + 2.*R2CSAxy;

    /* N Dipolar Interactions Contributions */
    R2CH = Dipolar_R2(omega_13C, omega_1H, w1, wr, pars->taus, pars->S2CHs, pars->tauf, pars->S2CHf, D_CH);
    R2CHr = Dipolar_R2(omega_13C, omega_1H, w1, wr, pars->taus, pars->S2CNs, pars->tauf, pars->S2CNf, D_CHr);
    R2CN = Dipolar_R2(omega_13C, omega_15N, w1, wr, pars->taus, pars->S2CNs, pars->tauf, pars->S2CNf, D_CN);
    R2CC = Dipolar_R2(omega_13C, omega_13C - wCOCa, w1, wr, pars->taus, pars->S2CCs, pars->tauf, pars->S2CCf, D_CC);
    return (Decimal) ((R2CSA + R2CH + R2CHr + R2CN + R2CC)*(Decimal)T_DOWN);
}
