//
// Created by ben on 16/03/2021.
//

#ifndef DYNAMIX_IMPACT_H
#define DYNAMIX_IMPACT_H

#include "../datatypes.h"

Decimal J0_IMPACT(Decimal omega, Decimal *A, Decimal tau_lb, Decimal tau_ub, unsigned int impact_n);
int TAU_IMPACT(Decimal *tau, Decimal tau_lb, Decimal tau_ub, unsigned int impact_n);
Decimal Dipolar_R1_IMPACT(Decimal omega_obs, Decimal omega_neigh, Decimal *A, Decimal tau_lb, Decimal tau_ub, unsigned int impact_n,
                          Decimal D);
Decimal Dipolar_R2_IMPACT(Decimal omega_obs, Decimal omega_neigh, Decimal w1, Decimal wr, Decimal *A, Decimal tau_lb, Decimal tau_ub, unsigned int impact_n, Decimal D);
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
);

Decimal Calc_15NR1_IMPACT(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m,
                          unsigned int mode);

Decimal Calc_15NR2_IMPACT(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m,
                          unsigned int mode);

Decimal Calc_13CR1_IMPACT(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m,
                          unsigned int mode);

Decimal Calc_13CR2_IMPACT(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m,
                          unsigned int mode);

#endif //DYNAMIX_MODEL_H
