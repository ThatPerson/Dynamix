//
// Created by ben on 16/03/2021.
//

#ifndef DYNAMIX_MODEL_H
#define DYNAMIX_MODEL_H

#include "../datatypes.h"

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
    Decimal Eas, Eaf;
    Decimal S2CHs, S2CHf;
    Decimal S2CCAcs, S2CCAcf;
    Decimal S2CCAps, S2CCApf;
    Decimal S2NHrs, S2NHrf;
    Decimal S2CHrs, S2CHrf;
    Decimal S2uf;
    Decimal papbS2, kex;
    Decimal dS2s, dS2f;
    Decimal GS2, Gtaur;

};

Decimal J0(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal S2uf);

Decimal J0_CC(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal S2uf);
Decimal PJ0(Decimal omega, Decimal r, Decimal ns, Decimal Gtaur);
//Decimal PJ0(Decimal omega, Decimal GS2, Decimal Gtaur);
Decimal Paramagnetic_R2(Decimal omega_N, Decimal omega_E, Decimal GS2, Decimal Conc, Decimal Gtaur, Decimal D, Decimal w1, Decimal wr);
Decimal Paramagnetic_R1(Decimal omega_N, Decimal omega_E, Decimal GS2, Decimal Conc, Decimal Gtaur, Decimal D);
/* R is effective distance, Conc is concentration */

Decimal Dipolar_R1(Decimal omega_obs, Decimal omega_neigh, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal S2uf,
           Decimal D);

Decimal Dipolar_R2(Decimal omega_obs, Decimal omega_neigh, Decimal w1, Decimal wr, Decimal taus, Decimal S2s, Decimal tauf,
           Decimal S2f, Decimal S2uf, Decimal D);

Decimal CSA_R2(Decimal omega, \
                         Decimal w1, \
                         Decimal wr, \
                         Decimal taus, \
                         Decimal S2s, \
                         Decimal tauf, \
                         Decimal S2f, \
                         Decimal S2uf, \
                         Decimal D2, \
                         Decimal  (*J_SD)(\
                            Decimal, \
                            Decimal, \
                            Decimal, \
                            Decimal, \
                            Decimal, \
                            Decimal)\
);

Decimal Calc_15NR1(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m, unsigned int mode);

Decimal Calc_15NR2(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m, unsigned int mode);

Decimal Calc_13CR1(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m, unsigned int mode);

Decimal Calc_13CR2(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m, unsigned int mode);

int AIMF_S2(Decimal order_params[3], struct Orient **A, Decimal *S2[], int length);

int GAF_S2(Decimal sig[3], struct Orient **A, struct Orient **B, Decimal *S2[], int length, unsigned int mode);

#endif //DYNAMIX_MODEL_H
