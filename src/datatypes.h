#ifndef DATATYPES_H
#define DATATYPES_H

typedef double Decimal;
typedef double _Complex Complex;

Decimal Dwig[5][5];
struct Model;
struct Orient;
struct Residue;
struct Relaxation;
struct rrargs;


void calculate_Y2(struct Orient * or);
void initialise_dwig(Decimal angle, Decimal Dw[5][5]);
void free_all(struct Model *m);
Decimal sq(Decimal x);
int sq_i(int x);
void rotate_Y2(struct Orient * or, Decimal alpha, Decimal beta, Decimal gamma);


#include "datatypes.c"
#endif
