#ifndef CROSEN_H
#define CROSEN_H

Decimal
simplex(Decimal (*func)(Decimal[], struct Residue *, struct Model *, unsigned int), Decimal start[], Decimal EPSILON,
        Decimal scale, struct Residue *resid, struct Model *mod);

#endif
