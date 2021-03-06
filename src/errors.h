#ifndef ERRORS_H
#define ERRORS_H

Decimal uniform_rand(void);
Decimal norm_rand(Decimal mean, Decimal std);
void calc_statistics(Decimal * vals, int length, Decimal * mean, Decimal * std);
int calc_errors(struct Model *m, unsigned int residue);


#include "errors.c"
#endif
