#ifndef CORRELATION_H
#define CORRELATION_H

int read_data_file(char *filename, struct Model * m);
int write_correlation_function_emf(char * fn, Decimal T, Decimal dT, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f);
int write_correlation_function_smf(char * fn, Decimal T, Decimal dT, Decimal tau, Decimal S2);
int main(int argc, char * argv[]);

#include "correlation.c"
#endif
