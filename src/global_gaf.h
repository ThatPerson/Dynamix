#ifndef GLOBAL_H
#define GLOBAL_H

Decimal optimize_global_chisq(Decimal *opts, struct Residue *r, struct Model *m, unsigned int params);

int run_global_iteration(struct Model *m, int i);

int calc_global_errors(struct Model *m);

int run_global(struct Model *m);

void global_fit_control(struct Model *m);

void global_errors_control(struct Model *m);

#endif
