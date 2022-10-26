//
// Created by ben on 16/03/2021.
//

#include <stdlib.h>
#include <complex.h>
#include "../chisq.h"
#include <stdio.h>
#include "model.h"
#include "impact.h"

Decimal J0_IMPACT(Decimal omega, Decimal *A, Decimal *impact_tau, unsigned int impact_n){

    Decimal J = 0;
    unsigned int i;
    for (i = 0; i < impact_n; i++) {
        Decimal taui = impact_tau[i];
        J += (Decimal) A[i] * (Decimal) taui / (1 + ((Decimal) omega * (Decimal) omega * (Decimal) taui * (Decimal) taui));
    }
    return J;
}

int TAU_IMPACT(Decimal *tau, Decimal tau_lb, Decimal tau_ub, unsigned int impact_n){
    // tau_lb and tau_ub are the timescales in ns.
    Decimal l_tau_lb = log(tau_lb);
    Decimal l_tau_ub = log(tau_ub);

    Decimal J = 0;
    Decimal inc = (l_tau_ub - l_tau_lb) / (impact_n - 1);
    unsigned int i;
    for (i = 0; i < impact_n; i++) {
        Decimal taui = l_tau_lb + i * inc;
        tau[i] = exp(l_tau_lb + i * inc);
    }
    return 1;
}

Decimal Dipolar_R1_IMPACT(Decimal omega_obs, Decimal omega_neigh, Decimal *A, Decimal *impact_tau, unsigned int impact_n,
           Decimal D) {
    // for SMF and SMFT models, set S2s = 0 and taus = 0.
    Decimal q = (0.1) * sq(D) * (\
        (J0_IMPACT(omega_neigh - omega_obs, A, impact_tau, impact_n)) \
 + 3 * (J0_IMPACT(omega_obs, A, impact_tau, impact_n)) \
 + 6 * (J0_IMPACT(omega_neigh + omega_obs, A, impact_tau, impact_n)) \
);
    return (Decimal) q;
}


Decimal Dipolar_R2_IMPACT(Decimal omega_obs, Decimal omega_neigh, Decimal w1, Decimal wr, Decimal *A, Decimal *impact_tau, unsigned int impact_n, Decimal D) {
    // for SMF and SMFT models, set S2s = 0 and taus = 0.
    return (Decimal) (\
        (1 / 20.) * sq(D) * (\
            (2 / 3.) * J0_IMPACT(2 * M_PI * (w1 + 2 * wr), A, impact_tau, impact_n) + \
            (2 / 3.) * J0_IMPACT(2 * M_PI * (w1 - 2 * wr), A, impact_tau, impact_n) + \
            (4 / 3.) * J0_IMPACT(2 * M_PI * (w1 + wr), A, impact_tau, impact_n) + \
            (4 / 3.) * J0_IMPACT(2 * M_PI * (w1 - wr), A, impact_tau, impact_n) + \
            (3.) * J0_IMPACT(omega_obs, A, impact_tau, impact_n) + \
            (1.) * J0_IMPACT(omega_neigh - omega_obs, A, impact_tau, impact_n) + \
            (6.) * J0_IMPACT(omega_neigh, A, impact_tau, impact_n) + \
            (6.) * J0_IMPACT(omega_neigh + omega_obs, A, impact_tau, impact_n)\
)\
);
}

/**
 * Calculates the R2 contribution of CSA of nuclei with frequency omega.
 * D2 is the squared D22x/D22y/D22xy variable (dependent on CSA)
 * Kurbanov2011, eqs 18-20.
 */

//Decimal *A, Decimal tau_lb, Decimal tau_ub, unsigned int impact_n
Decimal CSA_R2_IMPACT(Decimal omega, \
                         Decimal w1, \
                         Decimal wr, \
                         Decimal *A, \
                         Decimal tau_lb, \
                         Decimal tau_ub, \
                         unsigned int impact_n, \
                         Decimal D2, \
                         Decimal  (*J_SD)(\
                            Decimal, \
                            Decimal *, \
                            Decimal, \
                            Decimal, \
                            unsigned int
)\
) {
    return (Decimal) (\
        (1 / 45.) * (Decimal) D2 * (\
            (2 / 3.) * J_SD(2 * M_PI * (w1 - 2 * wr), A, tau_lb, tau_ub, impact_n) + \
            (2 / 3.) * J_SD(2 * M_PI * (w1 + 2 * wr), A, tau_lb, tau_ub, impact_n) + \
            (4 / 3.) * J_SD(2 * M_PI * (w1 - wr), A, tau_lb, tau_ub, impact_n) + \
            (4 / 3.) * J_SD(2 * M_PI * (w1 + wr), A, tau_lb, tau_ub, impact_n) + \
            (3.) * J_SD(omega, A, tau_lb, tau_ub, impact_n)\
)\
);
}

Decimal Calc_15NR1_IMPACT(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m, unsigned int mode, unsigned int modeRD) {
    (void) mode; // R1 has no wr/w1 dependence in this model.

    /* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
    Decimal field = relax->field * 1000000; // conversion to Hz

    /* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
     * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
     * software division) I've predefined it. Bond length taken as 1.02 A */
    unsigned int model = m->model;
    if (model != MOD_IMPACT) {
        ERROR("Should not be in IMPACT function with non IMPACT model\n");
        exit(-1);
    }


    Decimal omega_1H = T_DOWN * 2 * M_PI * field;
    Decimal omega_15N = T_DOWN * 2 * M_PI * field / 9.869683408806043;
    Decimal omega_13C = T_DOWN * 2 * M_PI * field / 3.976489314034722;
    Decimal omega_E = omega_1H * 658;
    Decimal d2x, d2y, d2xy, d2tot;

    Decimal *csa;
    csa = res->csaN;

    /* N CSA relaxation contribution */
    Decimal R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1NH, R1NHr, R1CN, R1CaN, R1E = 0;
    Decimal J1;


    //int bcpars_tdep(struct BCParameters *pars, struct BCParameters *npars, Decimal tfacts, Decimal tfactf) {

    struct BCParameters npars = *pars;

    d2tot = (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
    d2tot += (-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
    d2tot *= sq(0.000001 * omega_15N * T_UP);
    R1CSA = (2 / 15.) * d2tot * J0_IMPACT(omega_15N, pars->impact_a, m->impact_tau, m->impact_n);

    /* N Dipolar Interactions Contributions */
    R1NH = Dipolar_R1_IMPACT(omega_15N, omega_1H, pars->impact_a, m->impact_tau, m->impact_n, D_NH);
    R1NHr = 0;
    if (relax->hydrogen == PROTONATED) R1NHr = Dipolar_R1_IMPACT(omega_15N, omega_1H, pars->impact_a, m->impact_tau, m->impact_n, D_NHr);
    R1CN = Dipolar_R1_IMPACT(omega_15N, omega_13C, pars->impact_a, m->impact_tau, m->impact_n, D_CN);
    R1CaN = Dipolar_R1_IMPACT(omega_15N, omega_13C, pars->impact_a, m->impact_tau, m->impact_n, D_NCA);

    if (m->gd_mod == GD_MOD || m->gd_mod == GD_MOD_FIXTAU)
        R1E = Paramagnetic_R1(omega_15N, omega_E, npars.Gr6norm, relax->Gd, npars.Gtau, D_NE, relax->field);

    Decimal Rate = 0;
    if ((modeRD & REL_PARAMAG) != 0)
        Rate += R1E;

    if ((modeRD & REL_DYNAMIC) != 0)
        Rate += R1CSA + R1NH + R1NHr + R1CN + R1CaN;

    return Rate * T_DOWN;
    //return (Decimal) (R1CSA + R1NH + R1NHr + R1CN + R1CaN + R1E) * T_DOWN;
}



Decimal Calc_15NR2_IMPACT(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m, unsigned int mode, unsigned int modeRD) {
    /* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
    Decimal field = relax->field * 1000000; // conversion to Hz
    Decimal wr = relax->wr, w1 = relax->w1;
    if ((mode & COMPENSATE_W1) != 0)
        w1 = relax->compensate_w1;
    if ((mode & COMPENSATE_WR) != 0)
        wr = relax->compensate_wr;

    /* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
     * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
     * software division) I've predefined it. Bond length taken as 1.02 A */
    unsigned int model = m->model;
    if (model != MOD_IMPACT) {
        ERROR("Should not be in IMPACT function with non IMPACT model\n");
        exit(-1);
    }

    Decimal omega_1H = 2 * M_PI * field;
    Decimal omega_13C, omega_15N, omega_E;
    Decimal d2x, d2y, d2xy, d2tot;

    omega_13C = 2 * M_PI * field / 3.976489314034722;
    omega_15N = 2 * M_PI * field / 9.869683408806043;
    omega_E = omega_1H * 658;


    Decimal *csa;
    csa = res->csaN;

    d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
    d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
    d2xy = (Decimal) sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

    /* N CSA relaxation contribution */
    Decimal R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2NH, R2NHr, R2CN, R2CaN, R2E = 0;

    w1 *= T_DOWN;
    wr *= T_DOWN;
    omega_1H *= T_DOWN;
    omega_13C *= T_DOWN;
    omega_15N *= T_DOWN;
    omega_E *= T_DOWN;

    Decimal J0sum = 0;


    J0sum += (2 / 3.) * J0_IMPACT(2 * M_PI * (w1 - 2 * wr), pars->impact_a, m->impact_tau, m->impact_n);
    J0sum += (2 / 3.) * J0_IMPACT(2 * M_PI * (w1 + 2 * wr), pars->impact_a, m->impact_tau, m->impact_n);
    J0sum += (4 / 3.) * J0_IMPACT(2 * M_PI * (w1 - wr), pars->impact_a, m->impact_tau, m->impact_n);
    J0sum += (4 / 3.) * J0_IMPACT(2 * M_PI * (w1 + wr), pars->impact_a, m->impact_tau, m->impact_n);

    d2tot = (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
    d2tot += (-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
    d2tot *= sq(0.000001 * omega_15N * T_UP);
    R2CSA = (1 / 45.) * d2tot * (J0sum + 3 * J0_IMPACT(omega_15N, pars->impact_a, m->impact_tau, m->impact_n));


    /* N Dipolar Interactions Contributions */
    R2NH = Dipolar_R2_IMPACT(omega_15N, omega_1H, w1, wr, pars->impact_a, m->impact_tau, m->impact_n, D_NH);
    R2NHr = 0;
    if (relax->hydrogen == PROTONATED) R2NHr = Dipolar_R2_IMPACT(omega_15N, omega_1H, w1, wr, pars->impact_a, m->impact_tau, m->impact_n, D_NHr);
    R2CN = Dipolar_R2_IMPACT(omega_15N, omega_13C, w1, wr, pars->impact_a, m->impact_tau, m->impact_n, D_CN);
    R2CaN = Dipolar_R2_IMPACT(omega_15N, omega_13C, w1, wr, pars->impact_a, m->impact_tau, m->impact_n, D_NCA);

    if (m->gd_mod == GD_MOD || m->gd_mod == GD_MOD_FIXTAU)
        R2E = Paramagnetic_R2(omega_15N, omega_E, pars->Gr6norm, relax->Gd, pars->Gtau, D_NE, relax->field, w1, wr);

    Decimal Rate = 0;

    if ((modeRD & REL_PARAMAG) != 0) {
        Rate += R2E;
    }

    if ((modeRD & REL_DYNAMIC) != 0)
        Rate += R2CSA + R2NH + R2NHr + R2CN + R2CaN;

    return Rate * T_DOWN;

   // return (R2CSA + R2NH + R2NHr + R2CN + R2CaN  + R2E) * T_DOWN;
}

Decimal Calc_13CR1_IMPACT(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m, unsigned int mode, unsigned int modeRD) {
    (void) mode;
    /* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
    Decimal field = relax->field * 1000000; // conversion to Hz

    /* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
     * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
     * software division) I've predefined it. Bond length taken as 1.02 A */
    unsigned int model = m->model;
    if (model != MOD_IMPACT) {
        ERROR("Should not be in IMPACT function with non IMPACT model\n");
        exit(-1);
    }

    Decimal omega_1H = 2 * M_PI * field;
    Decimal omega_13C, omega_15N, wCOCa, omega_E;
    Decimal d2x, d2y, d2xy, d2tot;

    omega_13C = 2 * M_PI * field / 3.976489314034722;
    omega_15N = 2 * M_PI * field / 9.869683408806043;
    omega_E = omega_1H * 658; // electron frequency
    wCOCa = 120 * omega_13C * 0.000001;

    omega_1H *= T_DOWN;
    omega_13C *= T_DOWN;
    omega_15N *= T_DOWN;
    omega_E *= T_DOWN;
    wCOCa *= T_DOWN;

    Decimal *csa;
    csa = res->csaC;

/* N CSA relaxation contribution */
    Decimal R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1CH, R1CHr, R1CN, R1CCAc, R1CCAp, R1E = 0;
    Decimal J1;

    d2tot = (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
    d2tot += (-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
    d2tot *= sq(0.000001 * omega_13C * T_UP);
    R1CSA = (2 / 15.) * d2tot * J0_IMPACT(omega_13C, pars->impact_a, m->impact_tau, m->impact_n);



    /* N Dipolar Interactions Contributions */
    R1CH = Dipolar_R1_IMPACT(omega_13C, omega_1H, pars->impact_a, m->impact_tau, m->impact_n, D_CH);
    R1CHr = 0;
    if (relax->hydrogen == PROTONATED) R1CHr = Dipolar_R1_IMPACT(omega_13C, omega_1H, pars->impact_a, m->impact_tau, m->impact_n, D_CHr);
    R1CN = Dipolar_R1_IMPACT(omega_13C, omega_15N, pars->impact_a, m->impact_tau, m->impact_n, D_CN);
    R1CCAp = Dipolar_R1_IMPACT(omega_13C, omega_13C - wCOCa, pars->impact_a, m->impact_tau, m->impact_n, D_CCAp);
    R1CCAc = Dipolar_R1_IMPACT(omega_13C, omega_13C - wCOCa, pars->impact_a, m->impact_tau, m->impact_n, D_CCAc);

    if (m->gd_mod == GD_MOD || m->gd_mod == GD_MOD_FIXTAU)
        R1E = Paramagnetic_R1(omega_13C, omega_E, pars->Gr6norm, relax->Gd, pars->Gtau, D_CE, relax->field);

    Decimal Rate = 0;
    if ((modeRD & REL_PARAMAG) != 0)
        Rate += R1E;

    if ((modeRD & REL_DYNAMIC) != 0)
        Rate += R1CSA + R1CH + R1CHr + R1CN + R1CCAp + R1CCAc;

    return Rate * T_DOWN;

   // return (Decimal) (R1CSA + R1CH + R1CHr + R1CN + R1CCAp + R1CCAc + R1E) * T_DOWN;
}

Decimal Calc_13CR2_IMPACT(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m, unsigned int mode, unsigned int modeRD) {
    /* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
    Decimal field = relax->field * 1000000; // conversion to Hz
    Decimal wr = relax->wr, w1 = relax->w1;
    if ((mode & COMPENSATE_W1) != 0)
        w1 = relax->compensate_w1;
    if ((mode & COMPENSATE_WR) != 0)
        wr = relax->compensate_wr;
    /* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
     * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
     * software division) I've predefined it. Bond length taken as 1.02 A */
    unsigned int model = m->model;
    if (model != MOD_IMPACT) {
        ERROR("Should not be in IMPACT function with non IMPACT model\n");
        exit(-1);
    }

    Decimal omega_1H = 2 * M_PI * field;
    Decimal omega_13C, omega_15N, wCOCa, omega_E;
    Decimal d2x, d2y, d2xy, d2tot;

    omega_13C = 2 * M_PI * field / 3.976489314034722;
    omega_15N = 2 * M_PI * field / 9.869683408806043;
    omega_E = omega_1H * 658;
    wCOCa = 120 * omega_13C * 0.000001;

    w1 *= T_DOWN;
    wr *= T_DOWN;
    omega_1H *= T_DOWN;
    omega_13C *= T_DOWN;
    omega_15N *= T_DOWN;
    omega_E *= T_DOWN;
    wCOCa *= T_DOWN;

    Decimal *csa;
    csa = res->csaC;

    /* CSA relaxation contribution */
    Decimal R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2CH, R2CHr, R2CN, R2CCAc, R2CCAp, R2E = 0;
    Decimal J0sum = 0;

    J0sum +=
            (2 / 3.) * J0_IMPACT(2 * M_PI * (w1 - 2 * wr), pars->impact_a, m->impact_tau, m->impact_n);
    J0sum +=
            (2 / 3.) * J0_IMPACT(2 * M_PI * (w1 + 2 * wr), pars->impact_a, m->impact_tau, m->impact_n);
    J0sum += (4 / 3.) * J0_IMPACT(2 * M_PI * (w1 - wr), pars->impact_a, m->impact_tau, m->impact_n);
    J0sum += (4 / 3.) * J0_IMPACT(2 * M_PI * (w1 + wr), pars->impact_a, m->impact_tau, m->impact_n);

    d2tot = (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
    d2tot += (-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
    d2tot *= sq(0.000001 * omega_13C * T_UP);
    R2CSA = (1 / 45.) * d2tot *
            (J0sum + 3 * J0_IMPACT(omega_13C, pars->impact_a, m->impact_tau, m->impact_n));


    /* N Dipolar Interactions Contributions */
    R2CH = Dipolar_R2_IMPACT(omega_13C, omega_1H, w1, wr, pars->impact_a, m->impact_tau, m->impact_n, D_CH);
    R2CHr = 0;
    if (relax->hydrogen == PROTONATED) R2CHr = Dipolar_R2_IMPACT(omega_13C, omega_1H, w1, wr, pars->impact_a, m->impact_tau, m->impact_n, D_CHr);
    R2CN = Dipolar_R2_IMPACT(omega_13C, omega_15N, w1, wr, pars->impact_a, m->impact_tau, m->impact_n, D_CN);
    R2CCAp = Dipolar_R2_IMPACT(omega_13C, omega_13C - wCOCa, w1, wr, pars->impact_a, m->impact_tau, m->impact_n,
                        D_CCAp);
    R2CCAc = Dipolar_R2_IMPACT(omega_13C, omega_13C - wCOCa, w1, wr, pars->impact_a, m->impact_tau, m->impact_n,
                        D_CCAc);

    if (m->gd_mod == GD_MOD || m->gd_mod == GD_MOD_FIXTAU)
        R2E = Paramagnetic_R2(omega_13C, omega_E, pars->Gr6norm, relax->Gd, pars->Gtau, D_CE, relax->field, w1, wr);


    Decimal Rate = 0;
    if ((modeRD & REL_PARAMAG) != 0)
        Rate += R2E;

    if ((modeRD & REL_DYNAMIC) != 0)
        Rate += R2CSA + R2CH + R2CHr + R2CN + R2CCAp + R2CCAc;

    return Rate * T_DOWN;

    //return (Decimal) ((R2CSA + R2CH + R2CHr + R2CN + R2CCAp + R2CCAc + R2E) * (Decimal) T_DOWN);
}
