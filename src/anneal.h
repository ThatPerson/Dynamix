//
// Created by ben on 10/03/2021.
//

#ifndef DYNAMIX_ANNEAL_H
#define DYNAMIX_ANNEAL_H

#define MODE_RANDOM_RESTART 2
#define MODE_RANDOM_START 1
#define MODE_OUTPUT_FILE 4

Decimal propose(Decimal v, Decimal wobble);

Decimal anneal(
        Decimal (*func)(Decimal[], struct Residue *, struct Model *, unsigned int),
        Decimal parms[],
        const Decimal min[],
        const Decimal max[],
        unsigned int n_pars,
        Decimal T,
        unsigned int n_iter,
        Decimal wobble,
        Decimal thermostat,
        Decimal restart_prob,
        char filename[255],
        unsigned int mode,
        struct Residue *r,
        struct Model *m
);

#endif //DYNAMIX_ANNEAL_H
