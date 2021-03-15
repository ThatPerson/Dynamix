#ifndef SMF_H
#define SMF_H

#undef MOD_SMF
#undef MOD_SMFT
#define MOD_SMF 1
#define MOD_SMFT 4

Decimal J0_SMF(Decimal omega, Decimal tau, Decimal S2);
Decimal SMF_Dipolar_R1(Decimal omega_X, Decimal omega_Y, Decimal d, Decimal tau, Decimal S2);
Decimal SMF_Dipolar_R2(Decimal omega_X, Decimal omega_Y, Decimal d, Decimal tau, Decimal S2, Decimal J0sum);

Decimal SMF_R1(struct Residue *res, struct Relaxation* relax, Decimal tau, Decimal S2, unsigned int mode);
Decimal SMF_R2(struct Residue *res, struct Relaxation* relax, Decimal tau, Decimal S2, unsigned int mode);

#endif
