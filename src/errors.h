#ifndef ERRORS_H
#define ERRORS_H

long double uniform_rand(void);
double norm_rand(double mean, double std);
void calc_statistics(long double * vals, int length, long double * mean, long double * std);
void * calc_errors(void *input);


#include "errors.c"
#endif
