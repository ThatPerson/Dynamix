#ifndef DATATYPES_H
#define DATATYPES_H

typedef double Decimal;
typedef double _Complex Complex;

Decimal Dwig[5][5];
struct Model;
struct Orient;
struct Residue;
struct Relaxation;


void calculate_Y2(struct Orient * or);
void initialise_dwig(Decimal angle, Decimal Dw[5][5]);
void free_all(struct Model *m);
Decimal sq(Decimal x);
int sq_i(int x);
void rotate_Y2(struct Orient * or, Decimal alpha, Decimal beta, Decimal gamma);
Decimal uniform_rand(void);
void gen_params(Decimal *minv, Decimal *maxv, Decimal *pars, unsigned int n_pars);
void setup_paramlims(struct Model *m, struct Residue *r, Decimal * minv, Decimal * maxv);

#include "datatypes.c"
#endif
