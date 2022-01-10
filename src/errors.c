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
 * errors.c
 *
 * Provides functions to calculate errors, and then determine the mean and standard
 * deviation of these. Uses Monte-Carlo method for error analysis..
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "datatypes.h"
#include "chisq.h"
#include "crosen.h"


/**
 * Uses Box-Muller method to generate a normally distributed random number.
 * @param mean
 *  Mean of normal distribution to select random variable from
 * @param std
 *  Standard deviation of normal distribution from which random variable selected
 * @return Decimal
 *  Returns normally distributed random number
 */
Decimal norm_rand(Decimal mean, Decimal std) {
    // Box-Muller method
    Decimal rnd1, rnd2;
    rnd1 = (Decimal) uniform_rand();
    rnd2 = (Decimal) uniform_rand();
    Decimal unadj = sqrt(-2 * log(rnd1)) * cos(2 * M_PI * rnd2);

    //printf("%lf, %lf, %lf\n", rnd1, rnd2, mean + std*unadj);
    return mean + std * unadj;
}

/**
 * Calculates mean and standard deviation of values contained in array, then puts these into the given pointers.
 * @param vals
 *  Pointer to array of values to take statistics of
 * @param length
 *  Length of vals
 * @param mean
 *  Pointer to Decimal to contain mean
 * @param std
 *  Pointer to Decimal to contain standard deviation
 */
void calc_statistics(Decimal *vals, int length, Decimal *mean, Decimal *std) {
    int k;
    *mean = 0;
    Decimal m2 = 0;
    for (k = 0; k < length; k++) {
        *mean += vals[k];
        m2 += pow(vals[k], 2);
    }
    *mean = *mean / ((Decimal) length);
    m2 = m2 / ((Decimal) length);
    Decimal mdiff = m2 - pow(*mean, 2.);
    if (mdiff < 0 && fabsl(mdiff) < 0.000001) {
        ERROR("mdiff = %lf, making negative.", mdiff);
        // then all have converged to the same point and the std is 0.
        // but because Decimaling points are ew, it will give 'nan' for sqrt(-0.0000)
        // so we just invert the sign. It's not pretty but ah well.
        mdiff = -mdiff;
    }
    *std = sqrt(mdiff);

}

/**
 * Calculates errors for a given residue.
 * Error calculation is done by varying the relaxation values according to a normal distribution with mean (R) and standard deviation (Rerror/2).
 * Then simplex optimization is performed, and the newly optimized parameters stored. Then statistics of these back optimized parameters are taken,
 * with the errors being the standard deviation of these.
 * @param input
 *  rrarg struct containing the residue under consideration.
 */
int calc_errors(struct Model *m, unsigned int residue) {
    struct Residue *resid = &(m->residues[residue]);


    char outputdir[255];
    strcpy(outputdir, m->outputdir);

    //Decimal optim = resid->min_val;
    //char outputdir[255];
    //strcpy(outputdir, ((struct rrargs*)input)->outputdir);
    FILE *errp;
    char filename[300];
    sprintf(filename, "%s/errors_%d.dat", outputdir, residue + 1);
    errp = fopen(filename, "w");
    if (errp == NULL) {
        printf("%s not found.\n", filename);
        return -1;
    }

    unsigned int l, k, params;
    params = m->params;
    //printf("%d\n", params);
    Decimal *opts;
    opts = (Decimal *) malloc(sizeof(Decimal) * params);
    resid->errors_mean = (Decimal *) malloc(sizeof(Decimal) * params);
    resid->errors_std = (Decimal *) malloc(sizeof(Decimal) * params);
    resid->error_params = (Decimal **) malloc(sizeof(Decimal *) * params);
    for (k = 0; k < params; k++) {
        resid->error_params[k] = (Decimal *) malloc(sizeof(Decimal) * m->n_error_iter);
    }
    //resid->min_val = MIN_VAL;*/
    //Decimal val = 0;
    int p = 0;
    int ignore = -1;
    Decimal temp_R;

    // make backups of order parameters.
    resid->S2NHb = resid->S2NH;
    resid->S2CHb = resid->S2CH;
    resid->S2CNb = resid->S2CN;
    resid->S2CCb = resid->S2CC;

    struct BCParameters pars;
    int otb = opts_to_bcpars(resid->parameters, &pars, m, resid, &ignore);
    if (otb != 0) {
        free(opts);
        return -1;
    }
    int retries = 0;
    for (l = 0; l < m->n_error_iter; l++) {
        if (resid->ignore == 1) {
            return -1;
        }

        /* As the simplex optimization looks in the relaxation pointer, we temporarily
         * store the relaxation pointer elsewhere while we perform the optimization
         */
        resid->temp_relaxation = (resid->relaxation);
        resid->relaxation = NULL;
        resid->relaxation = (struct Relaxation *) malloc(sizeof(struct Relaxation) * resid->lim_relaxation);

        //if (i == 2)
        //	printf("Residue %d iteration %d: %lf\n ", i, l, opts[0]);

        for (k = 0; k < resid->n_relaxation; k++) {
            resid->relaxation[k].field = resid->temp_relaxation[k].field;
            resid->relaxation[k].wr = resid->temp_relaxation[k].wr;
            resid->relaxation[k].w1 = resid->temp_relaxation[k].w1;
            resid->relaxation[k].type = resid->temp_relaxation[k].type;
            resid->relaxation[k].compensate_wr = resid->temp_relaxation[k].compensate_wr;
            resid->relaxation[k].compensate_w1 = resid->temp_relaxation[k].compensate_w1;
            resid->relaxation[k].compensate = resid->temp_relaxation[k].compensate;
            resid->relaxation[k].compensate_wr = resid->temp_relaxation[k].compensate_wr;
            resid->relaxation[k].hydrogen = resid->temp_relaxation[k].hydrogen;
            resid->relaxation[k].Gd = resid->temp_relaxation[k].Gd;

/*struct Relaxation {


    int hydrogen;
    Decimal Gd;
};*/
            resid->relaxation[k].T = (m->t_error == 0)
                    ? resid->temp_relaxation[k].T
                    : norm_rand(resid->temp_relaxation[k].T, m->t_error / 2); // t_error is 2 standard deviations.


            temp_R = back_calc(resid, &(resid->temp_relaxation[k]), m, &ignore, &pars);
            //if (i == 2 && k == 2)
            //	printf("Iteration %d, R = %lf\n", l, temp_R); // should be the same each time
            resid->relaxation[k].R = norm_rand(temp_R, (resid->temp_relaxation[k].Rerror / 2.));


            LOG("%d %d %d prior: %lf, backcalc: %lf, new: %lf, error: %lf", i, l, k, resid->temp_relaxation[k].R,
                temp_R, resid->relaxation[k].R, resid->temp_relaxation[k].Rerror / 2.);

            //resid->relaxation[k].R = norm_rand(resid->temp_relaxation[k].R, (resid->temp_relaxation[k].Rerror)/2.);
            //printf("%lf, %lf -> %lf \n", resid->temp_relaxation[k].R, resid->temp_relaxation[k].Rerror, resid->relaxation[k].R);
            /* Rerror is Monte-Carlo calculated 2 standard deviations; norm_rand takes 1 std */
            resid->relaxation[k].Rerror = resid->temp_relaxation[k].Rerror;
        }

        resid->S2NH = norm_rand(resid->S2NHb, (resid->S2NHe / 2.)); // divided by two to get 1 standard deviation
        resid->S2CH = norm_rand(resid->S2CHb, (resid->S2CHe / 2.));
        resid->S2CN = norm_rand(resid->S2CNb, (resid->S2CNe / 2.));
        resid->S2CC = norm_rand(resid->S2CCb, (resid->S2CCe / 2.));

        for (k = 0; k < params; k++) {
            opts[k] = resid->parameters[k];
        }
        Decimal min = simplex(optimize_chisq, opts, 1.0e-16, 1, resid, m);
        LOG("%d %d min = %lf; orig = %lf", i, l, min, resid->min_val)
        // The actual value is more or less irrelevant, we're just interested in the values now in 'opts'
        //resid->error_params[l] = (Decimal *) malloc (sizeof(Decimal) * params);

        /* Return our pointers to how they were before... */
        free(resid->relaxation);
        resid->relaxation = NULL;
        resid->relaxation = resid->temp_relaxation;
        resid->temp_relaxation = NULL;


        if (min > 1000 * resid->min_val && retries < 1000) {
            // did not converge.
            retries++;
            ERROR("Error iteration %d for residue %d gave min_val = %lf (base mod = %lf), indicating no convergence. Rerunning iteration.\n",
                  l + 1, residue, min, resid->min_val);
            l--;
            continue;
        }
        fprintf(errp, "%d\t%lf", l + 1, min);
        for (k = 0; k < params; k++) {
            resid->error_params[k][p] = opts[k];
            fprintf(errp, "\t%le", opts[k]);
        }
        fprintf(errp, "\n");
        p++;


    }

    // restore of order parameters.
    resid->S2NH = resid->S2NHb;
    resid->S2CH = resid->S2CHb;
    resid->S2CN = resid->S2CNb;
    resid->S2CC = resid->S2CCb;

    resid->error_calcs = p;
    fclose(errp);
    free(opts);
    return 1;
}
