#ifndef CROSEN_H
#define CROSEN_H

double simplex(double (*func)(long double[], struct Residue*, unsigned int, unsigned int), long double start[],unsigned int n, long double EPSILON, long double scale, struct Residue * resid, unsigned int model, unsigned int or_variation);

#include "crosen.c"
#endif
