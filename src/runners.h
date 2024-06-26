//
// Created by ben on 08/03/2021.
//

#ifndef DYNAMIX_RUNNERS_H
#define DYNAMIX_RUNNERS_H

int run_residue(struct Model *m, unsigned int residue);

int run_errors(struct Model *m);

int print_backcalcs(struct Model *m);

int print_gaf(struct Model *m);

int print_residues(struct Model *m);

int print_errors(struct Model *m);

int run_fitting(struct Model *m);

#endif //DYNAMIX_RUNNERS_H
