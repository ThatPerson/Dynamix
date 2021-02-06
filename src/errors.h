#ifndef ERRORS_H
#define ERRORS_H

long double uniform_rand(void);
double norm_rand(double mean, double std);
void calc_statistics(long double * vals, int length, long double * mean, long double * std);
int calc_errors(struct Model *m, int residue);


#include "errors.c"
#endif
