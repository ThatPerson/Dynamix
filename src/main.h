#ifndef MAIN_H
#define MAIN_H

int run_residue(struct Model *m, int residue);
int run_errors(struct Model *m);
int main(int argc, char * argv[]);
int print_backcalcs(struct Model *m);
int print_gaf(struct Model *m);
int print_residues(struct Model *m);
int print_errors(struct Model *m);

#endif
