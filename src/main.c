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
 * main.c
 *
 * Implements main function which calls functions to initialise model, run model, and output.
 */

#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "datatypes.h"
#include "read_data.h"

#include "global_gaf.h"
#include "runners.h"


#include "verification.h"


/**
 * Start function. Initialises parameters, loads files, spawns threads and outputs data.
 * @opts filename
 *  Takes input as path to file containing System definition (*dx file)
 * @opts -e
 *  Enables error mode
 */
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    clock_t begin = clock();

    int numprocs, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    //printf("I am %d/%d.\n", myid, numprocs);

    /* Initialisation */
    start_time = time(0);
    srand((unsigned int) time(NULL) + myid);
    initialise_dwig(HALF_PI, Dwig);


    char system_file[255] = "";
    unsigned int i;
    int err_mod = DISABLED;
    for (i = 1; i < (unsigned int) argc; i++) {
        //printf("%s\n", argv[i]);
        if (strcmp(argv[i], "-e") == 0)
            err_mod = ENABLED;
        else if (strcmp(argv[i], "verify") == 0)
            verify_all();
        else
            strcpy(system_file, argv[i]);
    }

    if (strcmp(system_file, "") == 0) {
        printf("Please provide system file.\n");
        exit(-1);
    }


    struct Model m;
    int ret = read_system_file(system_file, &m);
    if (ret == -1) {
        printf("Error found, crashing peacefully...\n");
        free_all(&m);
        MPI_Finalize();
        exit(-1);
    }
    m.myid = myid;
    m.numprocs = numprocs;
    //printf("%d params.\n", m.params);
    //exit(-1);

    omp_set_num_threads(m.nthreads);

    m.error_mode = err_mod;
    if (m.error_mode == ENABLED && m.n_error_iter == 0) {
        printf("Please provide number of error iterations\n");
        ret = -1;
    }

    if (m.model == MOD_UNDEFINED) {
        printf("Model undefined\n");
        ret = -1;
    }

    if (ret == -1) {
        free_all(&m);
        MPI_Finalize();
    }


    //printf("%d\n", ret);
    if (m.myid == 0) {
        char filename[300];
        sprintf(filename, "%s/model.txt", m.outputdir);
        print_system(&m, filename);
    }

    if (m.global == LOCAL) {
        determine_residues(m.n_residues, myid, numprocs, &(m.proc_start), &(m.proc_end));
        printf("%d/%d: running %d - %d\n", myid + 1, numprocs, m.proc_start, m.proc_end);
        run_fitting(&m);
        printf("Worker %d: Printing residues...\n", m.myid + 1);
        print_residues(&m);
        if (m.error_mode == ENABLED) {
            printf("Worker %d: Running errors...\n", m.myid + 1);
            run_errors(&m);
            printf("Worker %d: Printing errors...\n", m.myid + 1);
            print_errors(&m);
        }
        printf("Worker %d: Printing gaf...\n", m.myid + 1);
        print_gaf(&m);
        print_backcalcs(&m);
    } else if (m.global == GLOBAL) {
        if (m.myid == 0) {
            global_fit_control(&m);
            printf("Worker %d: Fit control finished.\n", m.myid + 1);
        } else {
            determine_residues(m.n_residues, myid - 1, numprocs - 1, &(m.proc_start), &(m.proc_end));
            run_global(&m);
            printf("Worker %d: Printing residues...\n", m.myid + 1);
            print_residues(&m);
        }
        if (m.error_mode == ENABLED) {
            if (m.myid == 0) {
                global_errors_control(&m);
                printf("Worker %d: Error control finished.\n", m.myid + 1);
            } else {
                printf("Worker %d: Running errors...\n", m.myid + 1);
                calc_global_errors(&m);
                printf("Worker %d: Printing errors...\n", m.myid + 1);
                print_errors(&m);
            }
        }
        if (m.myid != 0) {
            printf("Worker %d: Printing gaf...\n", m.myid + 1);
            print_gaf(&m);
            print_backcalcs(&m);
        }
    }


    clock_t end = clock();

    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

    printf("Worker %d Finished (%lf seconds elapsed)\n", m.myid + 1, time_spent);

    free_all(&m);
    MPI_Finalize();
    return 0;
}
