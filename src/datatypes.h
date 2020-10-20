#ifndef DATATYPES_H
#define DATATYPES_H

long double Dwig[5][5];
struct Model;
struct Orient;
struct Residue;
struct Relaxation;
struct rrargs;


void calculate_Y2(struct Orient * or);
void initialise_dwig(double angle, long double Dw[5][5]);
void free_all(struct Model *m);
long double sq(long double x);
int sq_i(int x);
void rotate_Y2(struct Orient * or, double alpha, double beta, double gamma);


#include "datatypes.c"
#endif
