#ifndef CROSEN_H
#define CROSEN_H

double simplex(double (*func)(long double[], struct Residue*, struct Model *, unsigned int), long double start[], long double EPSILON, long double scale, struct Residue * resid, struct Model * mod);

#include "crosen.c"
#endif
