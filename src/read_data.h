#ifndef READ_DATA_H
#define READ_DATA_H

int read_resid_data(struct Model *m, char *filename, int dt);
int read_pp(struct Model *m, char *filename, int orient);
int read_relaxation_data(struct Model *m, char *filename);
int read_system_file(char *filename, struct Model * m);
int print_system(struct Model *m, char *filename);

#include "read_data.c"
#endif
