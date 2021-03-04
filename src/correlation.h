#ifndef CORRELATION_H
#define CORRELATION_H

int read_data_file(char *filename, struct Model * m);
int write_correlation_function_emf(char * fn, double T, double dT, double taus, double S2s, double tauf, double S2f);
int write_correlation_function_smf(char * fn, double T, double dT, double tau, double S2);
int main(int argc, char * argv[]);

#include "correlation.c"
#endif
