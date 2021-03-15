/*
 * Copyright (c) 2021 Ben Tatman, University of Warwick
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the Sosftware
 * is furnished to do so, subject to the following conditions:
 *
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software
 * - Any academic work deriving from the Software shall cite [CITATION].
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUR OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * anneal.c
 *
 * Implements simulated annealing of chsiq function. Much more effective at locating
 * global minima than Nelder-Mead simplex algorithm.
 */

#include <stdio.h>
#include <stdlib.h>
#include "datatypes.h"
#include "anneal.h"

Decimal propose(Decimal v, Decimal wobble) {
    Decimal rn = (Decimal) ((rand() % 100) - 50.) / 50.; // random number from -1. to +1.
    return v * (1 + rn * wobble);

}



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
)
{
    Decimal *best = (Decimal *) malloc(n_pars * sizeof(Decimal));
    Decimal *curr = (Decimal *) malloc(n_pars * sizeof(Decimal));
    Decimal *prop = (Decimal *) malloc(n_pars * sizeof(Decimal));
    FILE *fp = NULL;
    int file_open = 0;
    if ((mode & MODE_OUTPUT_FILE) != 0) {
        fp = fopen(filename, "w");
        if (fp == NULL)
            printf("Failed to open file %s\n", filename);
        else
            file_open = 1;
    }
    unsigned int k;
    if ((mode & MODE_RANDOM_START) != 0)
        gen_params(min, max, parms, n_pars);

    for (k = 0; k < n_pars; k++) {
        best[k] = parms[k];
        curr[k] = parms[k];
    }
    unsigned int itr;
    Decimal temp = T;
    Decimal cost_prop, cost_curr = func(curr, r, m, n_pars), cost_best = cost_curr;
    Decimal rnv;
    Decimal urand = 0;
    int jump = 0;
    for (itr = 0; itr < n_iter; itr++) {
        for (k = 0; k < n_pars; k++) {
            prop[k] = propose(curr[k], wobble);


        }
        temp /= thermostat;
        //temp = T * (1.*(n_iter - itr) / (1.*n_iter));
        cost_prop = func(prop, r, m, n_pars);


        urand = uniform_rand();
        if (cost_prop < cost_curr) {
            cost_curr = cost_prop;
            for (k = 0; k < n_pars; k++) {
                curr[k] = prop[k];
            }
            if (cost_prop < cost_best) {
                cost_best = cost_prop;
                for (k = 0; k < n_pars; k++) {
                    best[k] = prop[k];
                }
            }
        } else if (temp > 1e-10) {

            rnv = exp(1.*(cost_curr - cost_prop) / temp);
            if (rnv > urand) {
                jump = 1;
                cost_curr = cost_prop;
                for (k = 0; k < n_pars; k++) {
                    curr[k] = prop[k];
                }
            }
        }

        if ((mode & MODE_RANDOM_RESTART) != 0) {
            if (restart_prob > uniform_rand()) {
                gen_params(min, max, curr, n_pars);

                cost_curr = func(curr, r, m, n_pars);
            }

        }
        if (((mode & MODE_OUTPUT_FILE) != 0) && file_open == 1)
            fprintf(fp, "%d, %lf, %lf, %lf, %lf, %d\n", itr, temp, cost_curr, cost_prop, cost_best, jump);
        jump = 0;
        //printf("%d :: %lf :: %lf (%lf), %lf (%lf), %lf (%lf)\n", itr, temp, cost_curr, curr[0], cost_prop, prop[0], cost_best, best[0]);
    }
    for (k = 0; k < n_pars; k++) {
        parms[k] = best[k];
    }
    if (file_open == 1)
        fclose(fp);
    free(best);
    free(curr);
    free(prop);
    return cost_best;
}