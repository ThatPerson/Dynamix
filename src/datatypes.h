#ifndef DATATYPES_H
#define DATATYPES_H

long double Dwig[5][5];
struct Model;
struct Orient;
struct Residue;
struct Relaxation;
struct rrargs;


void calculate_Y2(struct Orient * or);
void initialise_dwig(void);
void free_all(struct Model *m);


#include "datatypes.c"
#endif
