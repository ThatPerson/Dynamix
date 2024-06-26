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
 * global_gaf.c
 *
 * Implements global model fitting (eg whole protein motions).
 * In theory by separating out a protein into different files you could simulate
 * secondary structure dynamics.
 * In the utils/global directory there is a python file to generate orientation files.
 */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "datatypes.h"
#include "model.h"
#include "chisq.h"
#include "anneal.h"
#include "crosen.h"
#include "errors.h"

Decimal optimize_global_chisq(Decimal *opts, struct Residue *r, struct Model *m, unsigned int params) {
    /* opts is a pointer to an array containing;
     *
     *  for SMF, [tau, S2]
     */
    (void) r;
    unsigned int l, c = 0;
    Decimal chisq = 0;
    for (l = 0; l < m->n_residues; l++) {
        //	Decimal optimize_chisq(Decimal * opts, struct Residue * resid, struct Model * m, unsigned int params) {
        if (m->residues[l].ignore == 1)
            continue;
        chisq += optimize_chisq(opts, &(m->residues[l]), m, params);
        c++;
    }

    return (Decimal) (chisq / c);
    /* normalise to number of relaxation measurements - otherwise when using like 85 the chisq becomes huge which hinders convergence */
}

int run_global_iteration(struct Model *m, unsigned int i) {
    struct Residue *resid = NULL;
    FILE *fp;
    char filename[300];
    sprintf(filename, "%s/fitting%d.dat", m->outputdir, m->myid);
    fp = fopen(filename, "a");
    if (fp == NULL) {
        printf("%s not found.\n", filename);
        return -1;
    }

    //printf("RESIDUE %d\n", i+1);
    //printf("Number of relaxations: %d\n", resid->n_relaxation);
    unsigned int k;
    unsigned int params = m->params;
    //printf("%d\n", params);
    Decimal *opts;
    opts = (Decimal *) malloc(sizeof(Decimal) * params);
    Decimal *anneal_pars = (Decimal *) malloc(sizeof(Decimal) * params);
    Decimal *minv = (Decimal *) malloc(sizeof(Decimal) * params);
    Decimal *maxv = (Decimal *) malloc(sizeof(Decimal) * params);
    if (maxv == NULL || minv == NULL)
        printf("Malloc failed.\n");
    setup_paramlims(m, 0.7, minv, maxv);

    /*if (resid->ignore == 1) {
        free(minv);
        free(maxv);
        free(anneal_pars);
        free(opts);
        fprintf(fp, "%d, %lf, -1, -1\n", i + 1, 1000.);
        fclose(fp);
        return -1;
    }*/

    unsigned int q;
    Decimal val;
    anneal(optimize_global_chisq, anneal_pars, minv, maxv, m->params, m->anneal_temp, 1000, m->anneal_wobb,
           m->anneal_therm, m->anneal_restart, NULL, MODE_RANDOM_RESTART + MODE_RANDOM_START, resid, m);

    //printf("%le, %le\n", opts[0], opts[1]);
    for (q = 0; q < m->n_nm_iter; q++) {
        for (k = 0; k < m->params; k++) {
            opts[k] = anneal_pars[k];
        }
        val = simplex(optimize_global_chisq, opts, 1.0e-16, 1, resid, m);
        if (val >= 1000000 || val < 0) {
            val = -1;
            for (k = 0; k < m->params; k++) {
                opts[k] = -1;
            }
        }

        fprintf(fp, "%d\t%lf", i + 1, val);
        for (k = 0; k < params; k++) {
            fprintf(fp, "\t%le", opts[k]);
        }
        fprintf(fp, "\n");
        unsigned int j;
        for (j = 0; j < m->n_residues; j++) {
            if (val < m->residues[j].min_val && val != -1) {
                //printf("New lowest %lf\n", val);
                m->residues[j].min_val = val;
                for (k = 0; k < params; k++) {
                    m->residues[j].parameters[k] = opts[k];
                }
            }
        }
    }

    printf("Control finished.\n");
    free(minv);
    free(maxv);


    free(opts);
    fclose(fp);
    return 1;
}

void global_errors_control(struct Model *m) {
    //Decimal
    //Decimal ** error_params;
    //m->residues[i].error_params[k][p]
    Decimal **paramstore;
    unsigned int i, l, k;
    paramstore = (Decimal **) malloc(sizeof(Decimal *) * m->n_error_iter);
    for (i = 0; i < m->n_error_iter; i++) {
        paramstore[i] = (Decimal *) malloc(sizeof(Decimal) * m->params);
    }

    int curr = 0, buff;

    int *finished = (int *) malloc(sizeof(int) * m->numprocs);
    for (l = 1; l < m->numprocs; l++) {
        finished[l] = 0;
    }
    int iter_per_proc = (m->n_error_iter / m->numprocs) + 3;
    int co = 0;
    for (l = 0; l < m->n_error_iter; l++) {
        co = 0;
        for (k = 1; k < m->numprocs; k++) {
            co += finished[k];
            if (finished[k] == 1)
                continue;
            MPI_Recv(&buff, 1, MPI_INT, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (buff == -1) {
                finished[k] = 1;
                continue;
            } else {
                MPI_Recv(paramstore[buff], m->params, MPI_DECIMAL, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        if (co == m->numprocs - 1)
            break;
    }

    //MPI_Send(&l, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    //MPI_Send(opts, m->params, MPI_DECIMAL, 0, 1, MPI_COMM_WORLD);

    for (i = 0; i < m->n_error_iter; i++) {
        for (k = 1; k < m->numprocs; k++) {
            MPI_Send(paramstore[i], m->params, MPI_DECIMAL, k, 2, MPI_COMM_WORLD);
        }
        free(paramstore[i]);
    }
    free(paramstore);
    free(finished);
}

int calc_global_errors(struct Model *m) {
    if (m->myid == 0)
        ERROR("I shouldn't reach this.\n");
    FILE *errp;
    char filename[300];
    sprintf(filename, "%s/errors%d.dat", m->outputdir, m->myid);
    errp = fopen(filename, "a");
    if (errp == NULL) {
        printf("%s not found.\n", filename);
        return -1;
    }

    unsigned int l, k, i;
    struct Residue *resid;
    Decimal *opts;
    opts = (Decimal *) malloc(sizeof(Decimal) * m->params);


    int p = 0;
    Decimal temp_R;
    int ignore = -1;

    for (i = 0; i < m->n_residues; i++) {
        resid = &(m->residues[i]);
        resid->S2NHb = resid->S2NH;
        resid->S2CHb = resid->S2CH;
        resid->S2CNb = resid->S2CN;
        resid->S2CCb = resid->S2CC;
        resid->errors_mean = (Decimal *) malloc(sizeof(Decimal) * m->params);
        resid->errors_std = (Decimal *) malloc(sizeof(Decimal) * m->params);
        resid->error_params = (Decimal **) malloc(sizeof(Decimal *) * m->params);
        for (k = 0; k < m->params; k++) {
            resid->error_params[k] = (Decimal *) malloc(sizeof(Decimal) * m->n_error_iter);
        }
    }


    struct BCParameters pars;
    unsigned int start, end;
    determine_residues(m->n_error_iter, m->myid - 1, m->numprocs - 1, &start, &end);

    for (l = start; l < end; l++) {

        printf("Iteration %d (%d-%d).\n", l, start, end - 1);
        for (i = 0; i < m->n_residues; i++) {
            resid = &(m->residues[i]);
            resid->temp_relaxation = (resid->relaxation);
            resid->relaxation = NULL;
            resid->relaxation = (struct Relaxation *) malloc(sizeof(struct Relaxation) * resid->lim_relaxation);

            int otb = opts_to_bcpars(resid->parameters, &pars, m, resid, &ignore);
            if (otb != 0) {
                free(opts);
                return -1;
            }

            for (k = 0; k < resid->n_relaxation; k++) {
                resid->relaxation[k].field = resid->temp_relaxation[k].field;
                resid->relaxation[k].wr = resid->temp_relaxation[k].wr;
                resid->relaxation[k].w1 = resid->temp_relaxation[k].w1;
                resid->relaxation[k].type = resid->temp_relaxation[k].type;
                resid->relaxation[k].T = resid->temp_relaxation[k].T;

                temp_R = back_calc(resid, &(resid->temp_relaxation[k]), m, &ignore, &pars);
                resid->relaxation[k].R = norm_rand(temp_R, (resid->temp_relaxation[k].Rerror / 2.));
                LOG("%d %d %d prior: %lf, backcalc: %lf, new: %lf, error: %lf", i, l, k, resid->temp_relaxation[k].R,
                    temp_R, resid->relaxation[k].R, resid->temp_relaxation[k].Rerror / 2.);
                resid->relaxation[k].Rerror = resid->temp_relaxation[k].Rerror;
            }
            resid->S2NH = norm_rand(resid->S2NHb, (resid->S2NHe / 2.)); // divided by two to get 1 standard deviation
            resid->S2CH = norm_rand(resid->S2CHb, (resid->S2CHe / 2.));
            resid->S2CN = norm_rand(resid->S2CNb, (resid->S2CNe / 2.));
            resid->S2CC = norm_rand(resid->S2CCb, (resid->S2CCe / 2.));
        }
        for (k = 0; k < m->params; k++) {
            opts[k] = m->residues[0].parameters[k]; // doesn't matter which we take it from, all are the same.
        }

        resid = NULL;
        Decimal min = simplex(optimize_global_chisq, opts, 1.0e-16, 1, resid, m);
        LOG("%d %d min = %lf", i, l, min);
        fprintf(errp, "%d\t%lf", l + 1, min);
        for (k = 0; k < m->params; k++) {
            for (i = 0; i < m->n_residues; i++)
                m->residues[i].error_params[k][p] = opts[k];
            fprintf(errp, "\t%le", opts[k]);
        }
        fprintf(errp, "\n");
        MPI_Send(&l, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(opts, m->params, MPI_DECIMAL, 0, 1, MPI_COMM_WORLD);
        p++;
        for (i = 0; i < m->n_residues; i++) {
            resid = &(m->residues[i]);
            free(resid->relaxation);
            resid->relaxation = NULL;
            resid->relaxation = resid->temp_relaxation;
            resid->temp_relaxation = NULL;
        }
    }
    for (i = 0; i < m->n_residues; i++) {
        m->residues[i].S2NH = m->residues[i].S2NHb;
        m->residues[i].S2CH = m->residues[i].S2CHb;
        m->residues[i].S2CN = m->residues[i].S2CNb;
        m->residues[i].S2CC = m->residues[i].S2CCb;
        m->residues[i].error_calcs = m->n_error_iter;
    }
    int finished = -1;
    MPI_Send(&finished, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    // Now receive
    // i - residue, k - parameter, p - iteration
    Decimal *buffer = (Decimal *) malloc(sizeof(Decimal) * m->params);
    Decimal sum = 0;
    for (p = 0; p < m->n_error_iter; p++) {
        //printf("%d listening for %d\n", m->myid, p);
        MPI_Recv(buffer, m->params, MPI_DECIMAL, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("%d received %d\n", m->myid, p);
        for (i = 0; i < m->n_residues; i++) {
            for (k = 0; k < m->params; k++) {
                //printf("%d, %d, %d (%d/%d)\n", p, i, k, m->myid+1, m->numprocs);
                m->residues[i].error_params[k][p] = buffer[k];
                sum = sum + buffer[k];
            }
        }

    }
    printf("%d/%d CHECKSUM: %f\n", m->myid + 1, m->numprocs, sum);
    // MPI_Send(paramstore[i], m->params, MPI_DECIMAL, k, 2, MPI_COMM_WORLD);

    fclose(errp);
    free(opts);
    free(buffer);
    return 1;
}

void global_fit_control(struct Model *m) {
    double *scores, *params;
    scores = (Decimal *) malloc(sizeof(Decimal) * m->numprocs);
    params = (Decimal *) malloc(sizeof(Decimal) * m->params);
    if (scores == NULL || params == NULL) {
        ERROR("Allocation failed!\n");
        free_all(m);
        exit(-1);
    }

    unsigned int i, min_proc = 1000;
    Decimal min = 1e9;
    Decimal d;
    for (i = 1; i < m->numprocs; i++) {
        MPI_Recv(&d, 1, MPI_DECIMAL, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        scores[i] = d;
        if (scores[i] < min) {
            min = scores[i];
            min_proc = i;
        }
    }

    int w = 1, l = 0;
    for (i = 1; i < m->numprocs; i++) {
        if (i == min_proc)
            MPI_Send(&w, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
        else
            MPI_Send(&l, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
        MPI_Send(&min, 1, MPI_DECIMAL, i, 4, MPI_COMM_WORLD);
    }
    if (min_proc == 1000) {
        ERROR("Did not find optimal!\n");
        free(m);
        free(params);
        free(scores);
        exit(-1);
    }
    MPI_Recv(params, m->params, MPI_DECIMAL, min_proc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (i = 1; i < m->numprocs; i++) {
        MPI_Send(params, m->params, MPI_DECIMAL, i, 3, MPI_COMM_WORLD);
    }
    free(params);
    free(scores);
}

/**
 * Operates residue optimization. Generates random parameter guesses and passes these to the simplex function.
 * @param input
 *  Pointer to rrarg containing thread information
 */
int run_global(struct Model *m) {
    if (m->numprocs < 2) {
        ERROR("Global GAF fitting requires at least two MPI processes (one control, one worker.\n");
        free_all(m);
        exit(-1);
    }

    unsigned int start, end;
    unsigned int l, k;
    int prodigal = -1;
    for (l = 0; l < m->n_residues; l++) {
        m->residues[l].parameters = (Decimal *) malloc(sizeof(Decimal) * m->params);
        m->residues[l].min_val = MIN_VAL;
        if (prodigal == -1 && m->residues[l].ignore == 0)
            prodigal = (int) l;
    }

    // launch
    if (m->myid == 0) {
        ERROR("I should not be here.\n");
    }
    determine_residues(m->n_anneal_iter, m->myid - 1, m->numprocs - 1, &start, &end);
    printf("%d/%d -> %d-%d (%d)\n", m->myid + 1, m->numprocs, start, end, m->n_anneal_iter);
    for (l = start; l < end; l++) {
        printf("\tRunning iteration %d\n", l);
        run_global_iteration(m, l);
    }
    // send score to control.
    // m->residues[j].parameters[k]
    Decimal min_v = m->residues[prodigal].min_val, true_min;
    printf("%d/%d: My score was %f, sending to Control.\n", m->myid + 1, m->numprocs, min_v);
    MPI_Send(&min_v, 1, MPI_DECIMAL, 0, 0, MPI_COMM_WORLD);
    int chosen_one = 0;
    MPI_Recv(&chosen_one, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&true_min, 1, MPI_DECIMAL, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (chosen_one == 1) {
        printf("%d/%d chosen.\n", m->myid + 1, m->numprocs);
        MPI_Send(m->residues[prodigal].parameters, (int) m->params, MPI_DECIMAL, 0, 2, MPI_COMM_WORLD);
    }
    MPI_Recv(m->residues[prodigal].parameters, (int) m->params, MPI_DECIMAL, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (l = 0; l < m->n_residues; l++) {
        if (l == prodigal)
            continue;
        for (k = 0; k < m->params; k++) {
            m->residues[l].min_val = true_min;
            m->residues[l].parameters[k] = m->residues[prodigal].parameters[k];
        }
    }
    return 1;
}

