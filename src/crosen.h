#ifndef CROSEN_H
#define CROSEN_H

double simplex(double (*func)(double[], struct Residue*, struct Model *, unsigned int), double start[], double EPSILON, double scale, struct Residue * resid, struct Model * mod);

#include "crosen.c"
#endif
