#ifndef ERRORS_H
#define ERRORS_H

double uniform_rand(void);
double norm_rand(double mean, double std);
void calc_statistics(double * vals, int length, double * mean, double * std);
int calc_errors(struct Model *m, int residue);


#include "errors.c"
#endif
