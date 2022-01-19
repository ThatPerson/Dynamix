//
// Created by ben on 19/01/2022.
//

#include "predict_rate.h"

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
#include "chisq.h"

int read_fit_data(struct Model *m) {
    FILE * fp = NULL;
    char filename[266] = "";
    sprintf(filename, "%s/%s", m->outputdir, (m->error_mode == ENABLED)?"errors.dat":"final.dat");
    fp = fopen(filename, "r");
    if (fp == NULL) {
        ERROR("Reading file failed");
        return -1;
    }
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *token;
    unsigned int curr = 0;
    unsigned int res = 0;
    unsigned int i;
    while ((read = getline(&line, &len, fp)) != -1) {
        // line is in line
        curr = 0;
        while ((token = strsep(&line, "\t"))) {
            if (curr == 0) {
                res = atoi(token);
                m->residues[res - 1].parameters = (Decimal *) malloc(sizeof(Decimal) * m->params);
                m->residues[res - 1].errors_std = (Decimal *) malloc(sizeof(Decimal) * m->params);
            } else if (curr > 2) {
                if (m->error_mode == ENABLED) {
                    if (curr % 2 != 0) {
                        // it's a parameter, specifically parameter
                        // (curr - 1) / 2
                        // curr = 3 should give 0,
                        // curr = 5 should give 1
                        // so ((curr - 3) / 2)
                        i = (curr - 3) / 2;
                        m->residues[res - 1].parameters[i] = atof(token);
                    } else {
                        // it's an error
                        i = (curr - 4) / 2;
                        m->residues[res - 1].errors_std[i] = atof(token) / 2;
                    }
                } else {
                    // it's a parameter
                    // curr - 3
                    // curr = 3 should give 0
                    // curr = 4 should give 1
                    i = curr - 3;
                    if (i > m->params)
                        continue;
                    m->residues[res - 1].parameters[i] = atof(token);
                    m->residues[res - 1].errors_std[i] = 0;

                }
            }
            curr++;
        }
    }
    fclose(fp);
    if (line) free(line);
    return 0;
}

void clean_relaxation(struct Model *m) {
    unsigned int i;
    for (i = 0; i < m->n_residues; i++) {
        free(m->residues[i].relaxation);
        m->residues[i].relaxation = (struct Relaxation *) malloc(sizeof(struct Relaxation) * N_RELAXATION);
        m->residues[i].lim_relaxation = N_RELAXATION;
        m->residues[i].n_relaxation = 0;
    }
}

int read_relaxation_set(struct Model *m, char *fn) {
    FILE * fp = NULL;
    fp = fopen(fn, "r");
    if (fp == NULL) {
        ERROR("Reading file failed");
        return -1;
    }
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, fp)) != -1) {
        line[strcspn(line, "\n")] = 0;
        read_relaxation_data(m, line);
    }
    fclose(fp);
    return 1;
}

int print_prate(struct Model *m) {
    char filename[300];
    unsigned int l;
    for (l = m->proc_start; l < m->proc_end; l++) {
        sprintf(filename, "%s/prate_%d.dat", m->outputdir, l + 1);
        back_calculate(&(m->residues[l]), m, filename, m->params);
    }
    return 1;
}

int main(int argc, char * argv[]) {
    // no point in bothering with MPI for this
    start_time = time(0);
    unsigned int seed = time(NULL);
    srand(seed);
    initialise_dwig(HALF_PI, Dwig);
    char system_file[255] = "";
    int err_mod = DISABLED;
    char relaxation_file[255] = "";
    unsigned int i;
    for (i = 1; i < (unsigned int) argc; i++) {
        if (strcmp(argv[i], "-e") == 0) {
            err_mod = ENABLED;
        } else if (strcmp(argv[i], "-p") == 0) {
            strcpy(relaxation_file, argv[i + 1]);
            i++;
        } else {
            strcpy(system_file, argv[i]);
        }
    }


    printf("=== Predict Rate ===\n");
    printf("Reading system parameters from %s\n", system_file);
    if (strcmp(system_file, "") == 0) {
        ERROR("No system file provided\n");
        exit(-1);
    }
    struct Model m;
    if (read_system_file(system_file, &m) == -1) {
        ERROR("Could not read system file\n");
        goto quit;
    }
    omp_set_num_threads(m.nthreads);
    determine_residues(m.n_residues, 0, 1, &(m.proc_start), &(m.proc_end));
    m.error_mode = err_mod;
    if (m.error_mode == ENABLED && m.n_error_iter == 0) {
        printf("Please provide number of error iterations\n");
        goto quit;
    }

    if (m.model == MOD_UNDEFINED) {
        printf("Model undefined\n");
        goto quit;
    }

    read_fit_data(&m);
    clean_relaxation(&m); // remove all previous relaxation data
    if (strcmp(relaxation_file, "") == 0) {
        ERROR("No relaxation file defined");
        goto quit;
    }

    if (read_relaxation_set(&m, relaxation_file) == -1) goto quit;

    char filename[266];
    sprintf(filename, "%s/prate.txt", m.outputdir);
    print_system(&m, filename);

    print_prate(&m);

    quit:
    free_all(&m);
    return 0;
}