/*
 * Copyright (c) 2021 Ben Tatman, University of Warwick
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the Sosftware
 * is furnished to do so, subject to the following conditions:
 *
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software
 * - Any academic work deriving from the Software shall cite [CITATION].
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUR OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * chisq.c
 *
 * Provides functions to backcalculate relaxation rate parameters, and calculate a
 * reduced chisq value from these to allow for parameter optimization.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "datatypes.h"
#include "models/model.h"
#include "errors.h"

int bcpars_init(struct BCParameters *pars, Decimal slow, Decimal fast, struct Model *m) {
    pars->taus = 0;
    pars->tauf = 0;
    pars->Eas = -1;
    pars->Eaf = -1;
    pars->S2NHs = slow;
    pars->S2NHf = fast;
    pars->S2NCSAxs = slow;
    pars->S2NCSAxf = fast;
    pars->S2NCSAys = slow;
    pars->S2NCSAyf = fast;
    pars->S2NCSAxys = slow;
    pars->S2NCSAxyf = fast;
    pars->S2CSAxs = slow;
    pars->S2CSAxf = fast;
    pars->S2CSAys = slow;
    pars->S2CSAyf = fast;
    pars->S2CSAxys = slow;
    pars->S2CSAxyf = fast;
    pars->S2CNs = slow;
    pars->S2CNf = fast;
    pars->S2CaNs = slow;
    pars->S2CaNf = fast;
    pars->S2CHs = slow;
    pars->S2CHf = fast;
    pars->S2CCAps = slow;
    pars->S2CCApf = fast;
    pars->S2CCAcs = slow;
    pars->S2CCAcf = fast;
    pars->S2NHrs = slow;
    pars->S2NHrf = fast;
    pars->S2CHrs = slow;
    pars->S2CHrf = fast;
    pars->S2uf = 1;
    pars->tau_uf = 0;
    pars->dS2s = 0;
    pars->dS2f = 0;
    pars->Gr6norm = 0;
    pars->Gtau = 0;
    unsigned int i;
    if (m->model == MOD_IMPACT) {
        for (i = 0; i < m->impact_n; i++)
            pars->impact_a[i] = 0;
    }

    return 0;
}

int bcpars_clean(struct BCParameters *pars, struct Model *m) {
    pars->S2NHs = 0; pars->S2NHf = 0;
    pars->S2NCSAxs = 0; pars->S2NCSAxf = 0;
    pars->S2NCSAys = 0; pars->S2NCSAyf = 0;
    pars->S2NCSAxys = 0; pars->S2NCSAxyf = 0;
    pars->S2CSAxs = 0; pars->S2CSAxf = 0;
    pars->S2CSAys = 0; pars->S2CSAyf = 0;
    pars->S2CSAxys = 0; pars->S2CSAxyf = 0;
    pars->S2CNs = 0; pars->S2CNf = 0;
    pars->S2CaNs = 0; pars->S2CaNf = 0;
    pars->taus = 0; pars->tauf = 0;
    pars->Eas = 0; pars->Eaf = 0;
    pars->S2CHs = 0; pars->S2CHf = 0;
    pars->S2CCAcs = 0; pars->S2CCAcf = 0;
    pars->S2CCAps = 0; pars->S2CCApf = 0;
    pars->S2NHrs = 0; pars->S2NHrf = 0;
    pars->S2CHrs = 0; pars->S2CHrf = 0;
    pars->S2uf = 0; pars->tau_uf = 0;
    pars->papbS2 = 0;
    pars->kex = 0;
    pars->dS2s = 0; pars->dS2f = 0;
    pars->Gr6norm = 0; pars->Gtau = 0;
    if (m->model == MOD_IMPACT) {
        pars->impact_a = (Decimal *) malloc(sizeof(Decimal) * m->impact_n);
        unsigned int i;
        for (i = 0; i < m->impact_n; i++) {
            pars->impact_a[i] = 0;
        }
    }

}

int bcpars_tdep(struct BCParameters *pars, struct BCParameters *npars, Decimal tfacts, Decimal tfactf) {
    npars->S2NHs = pars->S2NHs * tfacts;
    npars->S2NCSAxs = pars->S2NCSAxs * tfacts;
    npars->S2NCSAys = pars->S2NCSAys * tfacts;
    npars->S2NCSAxys = pars->S2NCSAxys * tfacts;
    npars->S2CSAxs = pars->S2CSAxs * tfacts;
    npars->S2CSAys = pars->S2CSAys * tfacts;
    npars->S2CSAxys = pars->S2CSAxys * tfacts;
    npars->S2CNs = pars->S2CNs * tfacts;
    npars->S2CaNs = pars->S2CaNs * tfacts;
    npars->S2CHs = pars->S2CHs * tfacts;
    npars->S2CCAps = pars->S2CCAps * tfacts;
    npars->S2CCAcs = pars->S2CCAcs * tfacts;
    npars->S2NHrs = pars->S2NHrs * tfacts;
    npars->S2CHrs = pars->S2CHrs * tfacts;

    npars->S2NHf = pars->S2NHf * tfactf;
    npars->S2NCSAxf = pars->S2NCSAxf * tfactf;
    npars->S2NCSAyf = pars->S2NCSAyf * tfactf;
    npars->S2NCSAxyf = pars->S2NCSAxyf * tfactf;
    npars->S2CSAxf = pars->S2CSAxf * tfactf;
    npars->S2CSAyf = pars->S2CSAyf * tfactf;
    npars->S2CSAxyf = pars->S2CSAxyf * tfactf;
    npars->S2CNf = pars->S2CNf * tfactf;
    npars->S2CaNf = pars->S2CaNf * tfactf;
    npars->S2CHf = pars->S2CHf * tfactf;
    npars->S2CCApf = pars->S2CCApf * tfactf;
    npars->S2CCAcf = pars->S2CCAcf * tfactf;
    npars->S2NHrf = pars->S2NHrf * tfactf;
    npars->S2CHrf = pars->S2CHrf * tfactf;
   
    return 0;
}

int bcpars_update(struct BCParameters *pars, Decimal slow, Decimal fast) {
    if (slow != -1) {
        pars->S2NHs = slow;
        pars->S2NCSAxs = slow;
        pars->S2NCSAys = slow;
        pars->S2NCSAxys = slow;
        pars->S2CSAxs = slow;
        pars->S2CSAys = slow;
        pars->S2CSAxys = slow;
        pars->S2CNs = slow;
        pars->S2CaNs = slow;
        pars->S2CHs = slow;
        pars->S2CCAps = slow;
        pars->S2CCAcs = slow;
        pars->S2NHrs = slow;
        pars->S2CHrs = slow;
    }
    if (fast != -1) {
        pars->S2NHf = fast;
        pars->S2NCSAxf = fast;
        pars->S2NCSAyf = fast;
        pars->S2NCSAxyf = fast;
        pars->S2CSAxf = fast;
        pars->S2CSAyf = fast;
        pars->S2CSAxyf = fast;
        pars->S2CNf = fast;
        pars->S2CaNf = fast;
        pars->S2CHf = fast;
        pars->S2CCApf = fast;
        pars->S2CCAcf = fast;
        pars->S2NHrf = fast;
        pars->S2CHrf = fast;
    }
    return 0;
}

void check_S2_violations(struct BCParameters *pars, int *violations, struct Model *m) {
    if (pars->S2NHs > 1 || pars->S2NHs < 0)
        (*violations)++;
    if (pars->S2NHf > 1 || pars->S2NHf < 0)
        (*violations)++;
    if (pars->S2CCAps > 1 || pars->S2CCAps < 0)
        (*violations)++;
    if (pars->S2CCApf > 1 || pars->S2CCApf < 0)
        (*violations)++;
    if (pars->S2CHs > 1 || pars->S2CHs < 0)
        (*violations)++;
    if (pars->S2CHf > 1 || pars->S2CHf < 0)
        (*violations)++;
    if (pars->S2CNs > 1 || pars->S2CNs < 0)
        (*violations)++;
    if (pars->S2CNf > 1 || pars->S2CNf < 0)
        (*violations)++;
    if (pars->S2uf > 1 || pars->S2uf < 0)
        (*violations)++;

    if (pars->papbS2 < 0)
        (*violations)++;

    //if (pars->Gtaur < pow(10, -7) || pars->Gtaur > pow(10, -5))
    if (pars->Gtau < 0) { (*violations)++; }
    if (pars->Gr6norm < 0.1) (*violations)++; // if Gr6norm < 0.1 then PJ0 is *10.

    if (m->model == MOD_IMPACT) {
        unsigned int i;
        for (i = 0; i < m->impact_n; i++) {
            if (pars->impact_a[i] < 0)
                (*violations)++;
        }
    }
    // For extreme values this can end up leading to overflow
}

int opts_to_bcpars(Decimal *opts, struct BCParameters *pars, struct Model *m, struct Residue *resid, int *violations) {
    Decimal S2s, S2f;
    int i;
    unsigned int model = m->model;
    Decimal *S2sP[] = {&(pars->S2NHs),
                       &(pars->S2NCSAxs),
                       &(pars->S2NCSAys),
                       &(pars->S2NCSAxys),
                       &(pars->S2CSAxs),
                       &(pars->S2CSAys),
                       &(pars->S2CSAxys),
                       &(pars->S2CNs),
                       &(pars->S2CaNs),
                       &(pars->S2CHs),
                       &(pars->S2CCAcs),
                       &(pars->S2CCAps)
    };
    Decimal *S2fP[] = {&(pars->S2NHf),
                       &(pars->S2NCSAxf),
                       &(pars->S2NCSAyf),
                       &(pars->S2NCSAxyf),
                       &(pars->S2CSAxf),
                       &(pars->S2CSAyf),
                       &(pars->S2CSAxyf),
                       &(pars->S2CNf),
                       &(pars->S2CaNf),
                       &(pars->S2CHf),
                       &(pars->S2CCAcf),
                       &(pars->S2CCApf)
    };
    struct Orient *As[] = {
            &(resid->orients[OR_NH]),
            &(resid->orients[OR_NCSAxx]),
            &(resid->orients[OR_NCSAyy]),
            &(resid->orients[OR_NCSAxx]),
            &(resid->orients[OR_CCSAxx]),
            &(resid->orients[OR_CCSAyy]),
            &(resid->orients[OR_CCSAxx]),
            &(resid->orients[OR_CN]),
            &(resid->orients[OR_NCA]),
            &(resid->orients[OR_CNH]),
            &(resid->orients[OR_CCAc]),
            &(resid->orients[OR_CCAp])
    };
    struct Orient *Bs[] = {
            &(resid->orients[OR_NH]),
            &(resid->orients[OR_NCSAxx]),
            &(resid->orients[OR_NCSAyy]),
            &(resid->orients[OR_NCSAyy]),
            &(resid->orients[OR_CCSAxx]),
            &(resid->orients[OR_CCSAyy]),
            &(resid->orients[OR_CCSAyy]),
            &(resid->orients[OR_CN]),
            &(resid->orients[OR_NCA]),
            &(resid->orients[OR_CNH]),
            &(resid->orients[OR_CCAc]),
            &(resid->orients[OR_CCAp])
    };
    bcpars_clean(pars, m);
    if (model == MOD_SMF || model == MOD_SMFT) {
        // tau, S2, [Ea]
        S2s = opts[1];
        S2f = (m->microsecond == ENABLED)?(resid->S2NH / S2s):1;
        bcpars_init(pars, S2s, S2f, m);
        pars->taus = opts[0];
        if (model == MOD_SMFT)
            pars->Eas = opts[2];

    } else if (model == MOD_EMF || model == MOD_DEMF || model == MOD_EMFT || model == MOD_DEMFT || model == MOD_RDEMFT || model == MOD_SDEMFT || model == MOD_SDEMF) {
        S2s = opts[1];
        S2f = resid->S2NH / S2s;
        if (model == MOD_DEMFT || model == MOD_DEMF || model == MOD_RDEMFT || model == MOD_SDEMFT || model == MOD_SDEMF) {
            S2f = opts[3];
        }
        bcpars_init(pars, S2s, S2f, m);

        pars->taus = opts[0];
        pars->tauf = opts[2];
        if (model == MOD_DEMFT || model == MOD_RDEMFT || model == MOD_SDEMFT) {
            pars->Eas = opts[4];
            pars->Eaf = opts[5];
        } else if (model == MOD_EMFT) {
            pars->Eas = opts[3];
            pars->Eaf = opts[4];
        }
        if (model == MOD_RDEMFT) {
            pars->papbS2 = opts[6];
            pars->kex = opts[7];
        } else if (model == MOD_SDEMFT) {
            pars->dS2s = opts[6];
            pars->dS2f = opts[7];
        } else if (model == MOD_SDEMF) {
            pars->dS2s = opts[4];
            pars->dS2f = opts[5];
        }
    } else if (model == MOD_GAF || model == MOD_GAFT) {
        bcpars_init(pars, 0, 0, m);
        pars->taus = opts[0];
        pars->tauf = opts[1];
        Decimal sigs[3] = {opts[2], opts[3], opts[4]};
        Decimal sigf[3] = {opts[5], opts[6], opts[7]};
        for (i = 0; i < 3; i++) {
            if (sigs[i] < 0 || sigs[i] > 0.52360)
                (*violations)++;
            if (sigf[i] < 0 || sigf[i] > 0.52360)
                (*violations)++;
        }
        GAF_S2(sigs, As, Bs, S2sP, 12, MODE_REAL);
        GAF_S2(sigf, As, Bs, S2fP, 12, MODE_REAL);
        if (model == MOD_GAFT) {
            pars->Eas = opts[8];
            pars->Eaf = opts[9];
        }
        pars->S2NHrs = 1;
        pars->S2NHrf = 1;
        pars->S2CHrs = 1;
        pars->S2CHrf = 1;
    } else if (model == MOD_BGF || model == MOD_BGFT) {
        bcpars_init(pars, 0, 0, m);
        pars->taus = opts[0];
        pars->tauf = opts[1];

        // taus, tauf, sigma_slow, alpha_slow, beta_slow, sigma_fast, alpha_fast, beta_fast

        Decimal sigs[3] = {0, 0, opts[2]};
        Decimal alphas = opts[3];
        Decimal betas = opts[4];
        Decimal sigf[3] = {0, 0, opts[5]};
        Decimal alphaf = opts[6];
        Decimal betaf = opts[7];
        for (i = 0; i < 3; i++) {
            if (sigs[i] < 0 || sigs[i] > 0.52360)
                (*violations)++;
            if (sigf[i] < 0 || sigf[i] > 0.52360)
                (*violations)++;
        }

        for (i = 0; i < N_OR; i++) {
            calculate_Y2(&(resid->orients[i]));
            rotate_Y2(&(resid->orients[i]), alphas, betas, 0);
        }
        GAF_S2(sigs, As, Bs, S2sP, 12, MODE_REAL);
        for (i = 0; i < N_OR; i++) {
            calculate_Y2(&(resid->orients[i]));
            rotate_Y2(&(resid->orients[i]), alphaf, betaf, 0);
        }
        GAF_S2(sigf, As, Bs, S2fP, 12, MODE_REAL);
        if (model == MOD_BGFT) {
            pars->Eas = opts[8];
            pars->Eaf = opts[9];
        }
        pars->S2NHrs = 1;
        pars->S2NHrf = 1;
        pars->S2CHrs = 1;
        pars->S2CHrf = 1;
    } else if (model == MOD_EGAF || model == MOD_EGAFT) {
        Decimal sigs[3] = {opts[2], opts[3], opts[4]};
        S2f = opts[5];
        for (i = 0; i < 3; i++) {
            if (sigs[i] < 0 || sigs[i] > 0.52360)
                (*violations)++;
        }
        bcpars_init(pars, 0, S2f, m);
        GAF_S2(sigs, As, Bs, S2sP, 12, MODE_REAL);
        pars->taus = opts[0];
        pars->tauf = opts[1];
        if (model == MOD_EGAFT) {
            pars->Eas = opts[6];
            pars->Eaf = opts[7];
        }
        pars->S2NHrs = 1;
        pars->S2NHrf = 1;
        pars->S2CHrs = 1;
        pars->S2CHrf = 1;
    } else if (model == MOD_BGAF || model == MOD_BGAFT) {

        Decimal sigs[3] = {opts[1], opts[2], opts[3]};
        bcpars_init(pars, 1, 1, m);
        GAF_S2(sigs, As, Bs, S2sP, 12, MODE_REAL);
        S2f = (m->microsecond == ENABLED) ? resid->S2NH / pars->S2NHs : 1;
        bcpars_update(pars, -1, S2f);
        pars->taus = opts[0];
        for (i = 0; i < 3; i++) {
            if (sigs[i] < 0 || sigs[i] > 0.52360)
                (*violations)++;
        }
        if (model == MOD_BGAFT) {
            pars->Eas = opts[4];
        }
        pars->S2NHrs = 1;
        pars->S2NHrf = 1;
        pars->S2CHrs = 1;
        pars->S2CHrf = 1;
    } else if (model == MOD_AIMF || model == MOD_AIMFT) {
        bcpars_init(pars, 0, 0, m);
        pars->taus = opts[0];
        pars->tauf = opts[1];
        Decimal sigs[3] = {opts[2], opts[3], opts[4]};
        Decimal sigf[3] = {opts[5], opts[6], opts[7]};
        AIMF_S2(sigs, As, S2sP, 12);
        AIMF_S2(sigf, As, S2fP, 12);
        if (model == MOD_AIMFT) {
            pars->Eas = opts[8];
            pars->Eaf = opts[9];
        }
        pars->S2NHrs = 1;
        pars->S2NHrf = 1;
        pars->S2CHrs = 1;
        pars->S2CHrf = 1;
    } else if (model == MOD_BAIMF || model == MOD_BAIMFT) {
        bcpars_init(pars, 1, 1, m);
        pars->taus = opts[0];
        Decimal sigs[3] = {opts[1], opts[2], opts[3]};
        AIMF_S2(sigs, As, S2sP, 12);
        S2f = (m->microsecond == ENABLED) ? resid->S2NH / pars->S2NHs : 1;
        bcpars_update(pars, -1, S2f);
        if (model == MOD_AIMFT) {
            pars->Eas = opts[4];
        }
        pars->S2NHrs = 1;
        pars->S2NHrf = 1;
        pars->S2CHrs = 1;
        pars->S2CHrf = 1;
    } else if (model == MOD_EAIMF || model == MOD_EAIMFT) {
        S2f = opts[5];
        bcpars_init(pars, 0, S2f, m);
        pars->taus = opts[0];
        pars->tauf = opts[1];
        Decimal sigs[3] = {opts[2], opts[3], opts[4]};
        AIMF_S2(sigs, As, S2sP, 12);
        if (model == MOD_EAIMFT) {
            pars->Eas = opts[6];
            pars->Eaf = opts[7];
        }
        pars->S2NHrs = 1;
        pars->S2NHrf = 1;
        pars->S2CHrs = 1;
        pars->S2CHrf = 1;
    } else if (m->model == MOD_IMPACT) {
        unsigned int n;
        for (n = 0; n < m->impact_n; n++) {
            pars->impact_a[n] = opts[n];
        }
    } else {
        ERROR("Model %d does not exist.\n", model);
        return 1;
    }
    if (m->ultrafast == ENABLED) {
        pars->S2uf = opts[m->UFS2];
        pars->tau_uf = opts[m->UFtau_uf];
    }
    if (m->gd_mod == ENABLED) {
        pars->Gr6norm = opts[m->GDS2];
        pars->Gtau = opts[m->GDtaur];
    }

//	Decimal upper_lim_tf = (Decimal) 1 * pow(10, -8 + T_S);
//	Decimal upper_lim_ts = (Decimal) 1 * pow(10, -5 + T_S);
//	Decimal lower_lim_tf = (Decimal) 1 * pow(10, -13 + T_S);

    Decimal upper_lim_tf = 1; // 10 ns
    Decimal upper_lim_ts = 10000; // 10 us
    Decimal lower_lim_tf = 0.0001; // 100 fs
    Decimal upper_lim_tuf = 1;

    if (m->microsecond == ENABLED) {
        upper_lim_tf = 1000; // 1 us
        upper_lim_ts = 100000; // 100 us
    }

    Decimal ttaus = pars->taus, ttauf = pars->tauf;
    if (pars->Eas != -1 || pars->Eaf != -1) {
        ttaus = temp_tau(pars->taus, pars->Eas, 300);
        ttauf = temp_tau(pars->tauf, pars->Eaf, 300);
    }
    check_S2_violations(pars, violations, m);
    // printf("Pars %lf (%lf)\n", ttaus, upper_lim_ts);
    if (ttaus < ttauf)
        (*violations)++;
    if (ttauf < pars->tau_uf)
        (*violations)++;
    if (ttaus < pars->tau_uf)
        (*violations)++;
    if (pars->tau_uf > upper_lim_tuf)
        (*violations)++;
    if (ttaus > upper_lim_ts)
        (*violations)++;
    if (ttauf > upper_lim_tf)
        (*violations)++;
    if (ttauf < lower_lim_tf && ttauf != 0) // ttauf = 1 for models with no fast timescale.
        (*violations)++;
    for (i = 0; i < 12; i++) {
        if (*(S2sP[0]) > 1 || (*S2sP[0]) < 0)
            (*violations)++;
        if (*(S2fP[0]) > 1 || (*S2fP[0]) < 0)
            (*violations)++;
    }
    if (pars->S2uf > 1 || pars->S2uf < 0.5)
        (*violations)++;
    return 0;
}


/**
 * Relaxation back calculation. Takes parameters, residue and relaxation and returns
 * calculated relaxation rate.
 * @params *opts
 *  Array of parameters
 * @params *resid
 *  Pointer to residue
 * @params *relax
 *  Pointer to relaxation rate
 * @params model
 *  Model definition
 * @params *violations
 *  Int returning number of violations of constraints (for chisq)
 * @return calculated R value
 */
Decimal back_calc(struct Residue *resid, struct Relaxation *relax, struct Model *m, const int *violations,
                  struct BCParameters *pars) {
    Decimal calc_R = -1;
    (void) violations;
    (void) m;
    //if (relax->R <= 0)
    //    return -1;

    Decimal (*C_15NR1)(struct Residue *, struct Relaxation *, struct BCParameters *, struct Model *, unsigned int) = &Calc_15NR1;
    Decimal (*C_15NR2)(struct Residue *, struct Relaxation *, struct BCParameters *, struct Model *, unsigned int) = &Calc_15NR2;
    Decimal (*C_13CR1)(struct Residue *, struct Relaxation *, struct BCParameters *, struct Model *, unsigned int) = &Calc_13CR1;
    Decimal (*C_13CR2)(struct Residue *, struct Relaxation *, struct BCParameters *, struct Model *, unsigned int) = &Calc_13CR2;

    if (m->model == MOD_IMPACT) {
        C_15NR1 = &Calc_15NR1_IMPACT;
        C_15NR2 = &Calc_15NR2_IMPACT;
        C_13CR1 = &Calc_13CR1_IMPACT;
        C_13CR2 = &Calc_13CR2_IMPACT;
    }

    switch (relax->type) {
        case R_15NR1:
            calc_R = C_15NR1(resid, relax, pars, m, NONE);
            if (relax->compensate != NO_COMPENSATE)
                calc_R -= C_15NR1(resid, relax, pars, m, relax->compensate);
            break;
        case R_15NR1p:
            calc_R = C_15NR2(resid, relax, pars, m, NONE);
            if (relax->compensate != NO_COMPENSATE)
                calc_R -= C_15NR2(resid, relax, pars, m, relax->compensate);
            break;
        case R_13CR1:
            calc_R = C_13CR1(resid, relax, pars, m, NONE);
            if (relax->compensate != NO_COMPENSATE)
                calc_R -= C_13CR1(resid, relax, pars, m, relax->compensate);
            break;
        case R_13CR1p:
            calc_R = C_13CR2(resid, relax, pars, m, NONE);
            if (relax->compensate != NO_COMPENSATE)
                calc_R -= C_13CR2(resid, relax, pars, m, relax->compensate);
            break;
        default:
            ERROR("Relaxation type %d unknown\n", relax->type);
            break;
    }

    return calc_R;
}

/**
 * Optimization function. Loops over all relaxation measurements and predicts using models.\n
 * Calculates (real - calc)^2 / (error^2) over all relaxation measurements and returns.\n
 * @param opts
 *  Pointer to array containing parameters.
 * @param resid
 *  Pointer to residue being considered
 * @param model
 *  MOD_SMF etc.
 * @return Returns chisq value.
 */
Decimal optimize_chisq(Decimal *opts, struct Residue *resid, struct Model *m, unsigned int params) {
    /* opts is a pointer to an array containing;
     *
     *  for SMF, [tau, S2]
     */

    unsigned int or_variations = m->or_variation;
    Decimal chisq = 0;
    int violations = 0;
    unsigned int i;
    Decimal alpha, beta, gamma;
    if (or_variations == VARIANT_A) {
        alpha = (Decimal) opts[m->OValpha];
        beta = (Decimal) opts[m->OVbeta];
        gamma = (Decimal) opts[m->OVgamma];
        if (fabs(alpha) > OR_LIMIT || fabs(beta) > OR_LIMIT || fabs(gamma) > OR_LIMIT)
            violations++;
        for (i = 0; i < N_OR; i++) {
            calculate_Y2(&(resid->orients[i]));
            rotate_Y2(&(resid->orients[i]), alpha, beta, gamma);
        }
    }

    struct BCParameters pars;
    int k = opts_to_bcpars(opts, &pars, m, resid, &violations);
    if (k != 0)
        return -1;
    int mvio = 0;
#pragma omp parallel for reduction(+:chisq, violations, mvio)
    for (i = 0; i < resid->n_relaxation; i++) {
        Decimal calc_R, mult = 1.;
        int l_vio = 0;
        if (resid->relaxation[i].R == -1)
            continue;
        calc_R = back_calc(resid, &(resid->relaxation[i]), m, &l_vio, &pars);
        violations += l_vio;
        if (m->cn_ratio == CNRATIO_ON) {
            mult = 1;
            if (resid->relaxation[i].type == R_13CR1 && resid->relaxation[i].type == R_13CR1p)
                mult = resid->cn;
        }
        Decimal bchisq, add;
        bchisq = chisq;
        add = mult * ((pow(resid->relaxation[i].R - calc_R, 2.)) / pow(resid->relaxation[i].Rerror, 2.));
        chisq += add;


        if (chisq > pow(10, 50) || chisq < 0.0000001) {
            // printf("Violation tripped: %le\n", chisq);

            mvio++;
        }

    }
    //dp += resid->n_relaxation;
    if (mvio != 0) {
        return 1e10;

    }

    chisq += 100000 * violations;

    //if (violations > 0) printf("%d violations\n", violations);
    if (violations > 0) {
        LOG("Chisq optimization had %d violations.", violations);
    }

    if (chisq < 0.0000001) {
        return 1e10;
    }


    Decimal S2NH, S2CH, S2CN, S2CC;

    unsigned int div = resid->n_relaxation;

    if (m->model == MOD_IMPACT) {
        Decimal sum_a = 0;
        unsigned int n;
        for (n = 0; n < m->impact_n; n++) {
            sum_a += pars.impact_a[n];
        }
        free(pars.impact_a);
        chisq += m->impact_w * ((pow((1 - resid->S2NH) - sum_a, 2.)) / pow(resid->S2NHe, 2.));
        div += m->impact_w;
    } else {
        S2NH = pars.S2NHs * pars.S2NHf * pars.S2uf;
        S2CH = pars.S2CHs * pars.S2CHf * pars.S2uf;
        S2CN = pars.S2CNs * pars.S2CNf * pars.S2uf;
        S2CC = pars.S2CCAps * pars.S2CCApf * pars.S2uf;

        if (m->WS2NH != 0) chisq += m->WS2NH * ((pow(resid->S2NH - S2NH, 2.)) / pow(resid->S2NHe, 2.));
        if (m->WS2CH != 0) chisq += m->WS2CH * ((pow(resid->S2CH - S2CH, 2.)) / pow(resid->S2CHe, 2.));
        if (m->WS2CN != 0) chisq += m->WS2CN * ((pow(resid->S2CN - S2CN, 2.)) / pow(resid->S2CNe, 2.));
        if (m->WS2CC != 0) chisq += m->WS2CC * ((pow(resid->S2CC - S2CC, 2.)) / pow(resid->S2CCe, 2.));
        div += m->WS2CC + m->WS2CH + m->WS2CN + m->WS2NH;
    }



    if (chisq < 0.0000001) // indicates overflow.
        return 1e10;
    return (Decimal) (chisq / div);
    /* normalise to number of relaxation measurements - otherwise when using like 85 the chisq becomes huge which hinders convergence */
}

/**
 * Back calculation function. Loops over all relaxation measurements and predicts using models.\n
 * Outputs file containing [relaxation id, calculated, real, real error].\n
 * @param opts
 *  Pointer to array containing parameters.
 * @param resid
 *  Pointer to residue being considered
 * @param model
 *  MOD_SMF etc.
 * @param filename
 *  File to output calculations into
 * @return Returns 1 if successful, else -1.
 */
int back_calculate(struct Residue *resid, struct Model *m, char *filename, unsigned int params) {
    /* opts is a pointer to an array containing;
     *
     *  for SMF, [tau, S2]
     */
    //if (opts[0] == -1)
    //	return -1;
    if (m->n_bc_iter <= 0)
        return 0;
    unsigned int or_variations = m->or_variation;
    Decimal calc_R;
    FILE *fp;
    fp = fopen(filename, "w");
    if (fp == NULL) {
        ERROR("%s not found.", filename);
        return -1;
    }
    int violations = 0;
    unsigned int i;
    Decimal alpha, beta, gamma;
    //(m->residues[l].parameters),
    Decimal *pars = resid->parameters;
    Decimal *pars_err = NULL;
    if (m->error_mode == ENABLED) pars_err = resid->errors_std;

    unsigned int c_bc_iter;

    Decimal *opts = NULL;
    opts = (Decimal *) malloc(sizeof(Decimal) * params);

    Decimal **rate_temp = NULL;
    rate_temp = (Decimal **) malloc(sizeof(Decimal *) * resid->n_relaxation);
    if (rate_temp == NULL)
        exit(-1);
    for (i = 0; i < resid->n_relaxation; i++) {
        rate_temp[i] = NULL;
        rate_temp[i] = (Decimal *) malloc(sizeof(Decimal) * m->n_bc_iter);
        if (rate_temp[i] == NULL)
            exit(-1);
    }



    for (c_bc_iter = 0; c_bc_iter < m->n_bc_iter; c_bc_iter++) {
        for (i = 0; i < params; i++) {
            opts[i] = (m->error_mode == ENABLED) ? norm_rand(pars[i], pars_err[i]) : pars[i];
        }
        if (or_variations == VARIANT_A) {
            alpha = (Decimal) opts[m->OValpha];
            beta = (Decimal) opts[m->OVbeta];
            gamma = (Decimal) opts[m->OVgamma];
            if (fabs(alpha) > OR_LIMIT || fabs(beta) > OR_LIMIT || fabs(gamma) > OR_LIMIT) {
                violations++;
            }
            for (i = 0; i < N_OR; i++) {
                calculate_Y2(&(resid->orients[i]));
                rotate_Y2(&(resid->orients[i]), alpha, beta, gamma);
            }
        }
        struct BCParameters pars;
        int k = opts_to_bcpars(opts, &pars, m, resid, &violations);
        if (k != 0) {
            return -1;
        }
        for (i = 0; i < resid->n_relaxation; i++) {
            calc_R = back_calc(resid, &(resid->relaxation[i]), m, &violations, &pars);
            rate_temp[i][c_bc_iter] = calc_R;
            /*fprintf(fp, "%d\t%lf\t%lf\t%lf", i, (calc_R < 0 ? -1. : calc_R), resid->relaxation[i].R,
                    resid->relaxation[i].Rerror);

            if (VERBOSE)
                fprintf(fp, "\t%lf\t%lf\t%lf\t%lf\t%d", resid->relaxation[i].field, resid->relaxation[i].wr,
                        resid->relaxation[i].w1, resid->relaxation[i].T, resid->relaxation[i].type);
            fprintf(fp, "\n");*/
        }
    }

    free(opts);
    Decimal calcR, calcE;
    for (i = 0; i < resid->n_relaxation; i++) {
        calc_statistics(rate_temp[i], m->n_bc_iter, &calcR, &calcE);
        fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf", i, calcR, 2 * calcE, resid->relaxation[i].R,
                    resid->relaxation[i].Rerror);

        if (VERBOSE)
            fprintf(fp, "\t%lf\t%lf\t%lf\t%lf\t%d\t%f", resid->relaxation[i].field, resid->relaxation[i].wr,
                  resid->relaxation[i].w1, resid->relaxation[i].T, resid->relaxation[i].type, resid->relaxation[i].Gd);
        fprintf(fp, "\n");
        free(rate_temp[i]);
    }
    free(rate_temp);

    fclose(fp);
    return 1;
}
