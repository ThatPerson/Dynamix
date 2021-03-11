//
// Created by ben on 08/03/2021.
//

//#include "runners.h"



/**
 * Operates residue optimization. Generates random parameter guesses and passes these to the simplex function.
 * @param input
 *  Pointer to rrarg containing thread information
 */
int run_residue(struct Model *m, unsigned int residue) {
    struct Residue * resid = &(m->residues[residue]);



    unsigned int model = m->model;
    char outputdir[255];
    strcpy(outputdir, m->outputdir);
    FILE *fp;
    char filename[300];

    sprintf(filename, "%s/residue_%d.dat", outputdir, residue+1);
    fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("%s not found.\n", filename);
        return -1;
    }


    //printf("RESIDUE %d\n", i+1);
    //printf("Number of relaxations: %d\n", resid->n_relaxation);
    unsigned int l, k, q;
    unsigned int params = m->params;
    //printf("%d\n", params);
    Decimal *opts;
    opts = (Decimal *) malloc (sizeof(Decimal) * params);
    Decimal *anneal_pars = (Decimal *) malloc(sizeof(Decimal) * params);
    resid->parameters = (Decimal *) malloc (sizeof(Decimal) * params);
    if (resid->parameters == NULL || opts == NULL) {
        printf("Allocation failed.\n");
    }
    resid->min_val = MIN_VAL;

    Decimal *minv = (Decimal *) malloc(sizeof(Decimal) * params);
    Decimal *maxv = (Decimal *) malloc(sizeof(Decimal) * params);

    setup_paramlims(m, resid, minv, maxv);
    Decimal val;

    if (resid->ignore == 1) {
        free(minv);
        free(maxv);
        free(anneal_pars);
        free(opts);
        fprintf(fp, "%d, %lf, -1, -1\n", residue + 1, 1000.);
        fclose(fp);
        return -1;
    }

    for (l = 0; l < m->n_anneal_iter; l++) {
        val = anneal(optimize_chisq, anneal_pars, minv, maxv, params, m->anneal_temp, 1000, m->anneal_wobb, m->anneal_therm, m->anneal_restart, NULL, MODE_RANDOM_START +MODE_RANDOM_RESTART, resid, m);
        for (q = 0; q < m->n_nm_iter; q++) {
            for (k = 0; k < params; k++) {
                opts[k] = anneal_pars[k];
            }
            Decimal val = simplex(optimize_chisq, opts, 1.0e-16, 1, resid, m);
            if (val >= 1000000 || val < 0) {
                val = -1;
                for (k = 0; k < params; k++) {
                    opts[k] = -1;
                }
            }
            fprintf(fp, "%d\t%lf", residue + 1, val);
            for (k = 0; k < params; k++) {
                fprintf(fp, "\t%le", opts[k]);
            }

            if (val < resid->min_val && val != -1) {
                //printf("New lowest %lf\n", val);
                resid->min_val = val;
                for (k = 0; k < params; k++) {
                    resid->parameters[k] = opts[k];
                }
            }
            fprintf(fp, "\n");
        }
    }
    free(minv);
    free(maxv);
    free(anneal_pars);
    free(opts);
    fclose(fp);
    return 0;
}

int run_fitting(struct Model *m) {
    unsigned int i;
    for (i = m->proc_start; i < m->proc_end; i++) {
        printf("\tWorker %d/%d: Residue %d.\n", m->myid + 1, m->numprocs, i + 1);
        run_residue(m, i);
    }
    return 1;
}

int run_errors(struct Model *m) {
    unsigned int i;
    for (i = m->proc_start; i < m->proc_end; i++) {
        printf("\tWorker %d/%d: Residue %d (errors).\n", m->myid + 1, m->numprocs, i+1);
        calc_errors(m, i);
    }
    return 1;
}

int print_errors(struct Model *m) {
    FILE * ep = NULL;
    char filename[300];
    unsigned int l, k, i;
    sprintf(filename, "%s/errors_proc%d.dat", m->outputdir, m->myid);
    ep = fopen(filename, "w");
    if (ep == NULL) {
        ERROR("%s not found.", filename);
        free_all(m);
        return -1;
    }
    for (l = m->proc_start; l < m->proc_end; l++) {
        for (k = 0; k < m->params; k++) {
            calc_statistics(m->residues[l].error_params[k], m->residues[l].error_calcs, &(m->residues[l].errors_mean[k]), &(m->residues[l].errors_std[k]));
        }

        fprintf(ep, "%d\t%lf\t%lf", l+1, m->residues[l].S2NH, m->residues[l].min_val);

        for (i = 0; i < m->params; i++) {

            fprintf(ep, "\t%le\t%le", m->residues[l].parameters[i], 2 * m->residues[l].errors_std[i]);
            /* WARNING: I'm printing the actual minimized parameters with the errors from calculation.
             * The error_means are generally not minimal. */
        }

        fprintf(ep, "\n");
    }
    fclose(ep);
    return 1;
}

int print_residues(struct Model *m) {
    FILE * fp;
    char filename[300];
    unsigned int l, i;
    sprintf(filename, "%s/final_proc%d.dat", m->outputdir, m->myid);
    fp = fopen(filename, "w");
    if (fp == NULL) {
        ERROR("%s not found.", filename);
        free_all(m);
        return -1;
    }
    for (l = m->proc_start; l < m->proc_end; l++) {
        if (m->residues[l].min_val == MIN_VAL) {
            m->residues[l].min_val = -1.;
            for (i = 0; i < m->params; i++) {
                m->residues[l].parameters[i] = -1.;
            }
        }
        fprintf(fp, "%d\t%lf\t%lf", l+1, m->residues[l].S2NH, m->residues[l].min_val);
        for (i = 0; i < m->params; i++) {
            fprintf(fp, "\t%le", m->residues[l].parameters[i]);
        }
        fprintf(fp, "\t%d\n", m->residues[l].error_calcs);
    }
    fclose(fp);
    return 1;

}

int print_gaf(struct Model *m) {
    FILE * gaf;
    unsigned int l, i;
    char filename[300];
    sprintf(filename, "%s/gaf_proc%d.dat", m->outputdir, m->myid);
    gaf = fopen(filename, "w");
    if (gaf == NULL) {
        ERROR("%s not found.", filename);
        free_all(m);
        return -1;
    }
    FILE * orderparams;
    sprintf(filename, "%s/orderparams_proc%d.dat", m->outputdir, m->myid);
    orderparams = fopen(filename, "w");
    if (orderparams == NULL) {
        ERROR("%s not found.", filename);
        free_all(m);
        return -1;
    }

    for (l = m->proc_start; l < m->proc_end; l++) {
        Decimal alpha, beta, gamma;
        if (m->or_variation == VARIANT_A) {
            alpha = (Decimal) m->residues[l].parameters[m->params-3];
            beta = (Decimal) m->residues[l].parameters[m->params-2];
            gamma = (Decimal) m->residues[l].parameters[m->params-1];
            for (i = 0; i < N_OR; i++) {
                calculate_Y2(&(m->residues[l].orients[i]));
                rotate_Y2(&(m->residues[l].orients[i]), alpha, beta, gamma);
            }
        }

        /* I really need to clean this mess up at some point... */

        if ((m->model == MOD_GAF || m->model == MOD_GAFT) && gaf != NULL) {
            Decimal S2NHs, S2NHf, S2CHs, S2CHf, S2CNs, S2CNf, S2CCs, S2CCf;
            Decimal sigs[3] = {m->residues[l].parameters[2], m->residues[l].parameters[3], m->residues[l].parameters[4]};
            Decimal sigf[3] = {m->residues[l].parameters[5], m->residues[l].parameters[6], m->residues[l].parameters[7]};
            struct Orient *Os[] = {&(m->residues[l].orients[OR_NH]), &(m->residues[l].orients[OR_CNH]), &(m->residues[l].orients[OR_CN]), &(m->residues[l].orients[OR_CCAp])};
            Decimal *S2s[] = {&S2NHs, &S2CHs, &S2CNs, &S2CCs};
            Decimal *S2f[] = {&S2NHf, &S2CHf, &S2CNf, &S2CCf};

            GAF_S2(sigs, Os, Os, S2s, 4, MODE_REAL);
            GAF_S2(sigf, Os, Os, S2f, 4, MODE_REAL);
            fprintf(orderparams, "%d\t", l+1);
            fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2NH, m->residues[l].S2NHe);
            if (m->WS2CH != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CHs*S2CHf, m->residues[l].S2CH, m->residues[l].S2CHe);
            if (m->WS2CN != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CNs*S2CNf, m->residues[l].S2CN, m->residues[l].S2CNe);
            if (m->WS2CC != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CCs*S2CCf, m->residues[l].S2CC, m->residues[l].S2CCe);
            fprintf(orderparams, "\n");

            Decimal taus = m->residues[l].parameters[0];
            Decimal tauf = m->residues[l].parameters[1];
            if (m->model == MOD_GAFT) {
                Decimal Eas = m->residues[l].parameters[8];
                Decimal Eaf = m->residues[l].parameters[9];
                taus = temp_tau(taus, Eas, 300);
                tauf = temp_tau(tauf, Eaf, 300);
            }

            fprintf(gaf, "%d\t%le\t%lf\t%le\t%lf\n", l+1, taus, S2NHs, tauf, S2NHf);
        } else if ((m->model == MOD_AIMF || m->model == MOD_AIMFT) && gaf != NULL) {
            Decimal S2NHs, S2NHf, S2CHs, S2CHf, S2CNs, S2CNf, S2CCs, S2CCf;
            Decimal Ss[3] = {m->residues[l].parameters[2], m->residues[l].parameters[3], m->residues[l].parameters[4]};
            Decimal Sf[3] = {m->residues[l].parameters[5], m->residues[l].parameters[6], m->residues[l].parameters[7]};
            struct Orient *Os[] = {&(m->residues[l].orients[OR_NH]), &(m->residues[l].orients[OR_CNH]), &(m->residues[l].orients[OR_CN]), &(m->residues[l].orients[OR_CCAp])};
            Decimal *S2s[] = {&S2NHs, &S2CHs, &S2CNs, &S2CCs};
            Decimal *S2f[] = {&S2NHf, &S2CHf, &S2CNf, &S2CCf};

            AIMF_S2(Ss, Os, S2s, 4);
            AIMF_S2(Sf, Os, S2f, 4);
            fprintf(orderparams, "%d\t", l+1);
            fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2NH, m->residues[l].S2NHe);
            if (m->WS2CH != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CHs*S2CHf, m->residues[l].S2CH, m->residues[l].S2CHe);
            if (m->WS2CN != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CNs*S2CNf, m->residues[l].S2CN, m->residues[l].S2CNe);
            if (m->WS2CC != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2CCs*S2CCf, m->residues[l].S2CC, m->residues[l].S2CCe);
            fprintf(orderparams, "\n");

            Decimal taus = m->residues[l].parameters[0];
            Decimal tauf = m->residues[l].parameters[1];
            if (m->model == MOD_AIMFT) {
                Decimal Eas = m->residues[l].parameters[8];
                Decimal Eaf = m->residues[l].parameters[9];
                taus = temp_tau(taus, Eas, 300);
                tauf = temp_tau(tauf, Eaf, 300);
            }

            fprintf(gaf, "%d\t%le\t%lf\t%le\t%lf\n", l+1, taus, S2NHs, tauf, S2NHf);
        } else if ((m->model == MOD_EGAF || m->model == MOD_EGAFT) && gaf != NULL) {
            Decimal S2NHs, S2NHf = (Decimal) m->residues[l].parameters[5];
            Decimal *S2[] = {&S2NHs};
            /* Approximate as just the S2NHs and S2NHf */
            struct Orient *As[] = {&(m->residues[l].orients[OR_NH])};
            Decimal sigs[3] = {m->residues[l].parameters[2], m->residues[l].parameters[3], m->residues[l].parameters[4]};
            GAF_S2(sigs, As, As, S2, 1, MODE_REAL);


            fprintf(orderparams, "%d\t", l+1);
            fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2NH, m->residues[l].S2NHe);
            if (m->WS2CH != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2CH, m->residues[l].S2CHe);
            if (m->WS2CN != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2CN, m->residues[l].S2CNe);
            if (m->WS2CC != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2CC, m->residues[l].S2CCe);
            fprintf(orderparams, "\n");


            Decimal taus = m->residues[l].parameters[0];
            Decimal tauf = m->residues[l].parameters[1];
            if (m->model == MOD_EGAFT) {
                Decimal Eas = m->residues[l].parameters[6];
                Decimal Eaf = m->residues[l].parameters[7];
                taus = temp_tau(taus, Eas, 300);
                tauf = temp_tau(tauf, Eaf, 300);
            }
            fprintf(gaf, "%d\t%le\t%lf\t%le\t%lf\n", l+1, taus, S2NHs, tauf, S2NHf);
            /* This gives rise to an uninitialized warning, but this is initialized on line 445 inside another if so it should be fine. */
        } else if (m->model == MOD_DEMF || m->model == MOD_DEMFT) {
            Decimal S2NHs, S2NHf;
            S2NHs = (Decimal) m->residues[l].parameters[1];
            S2NHf = (Decimal) m->residues[l].parameters[3];
            fprintf(orderparams, "%d\t", l+1);
            fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2NH, m->residues[l].S2NHe);
            if (m->WS2CH != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2CH, m->residues[l].S2CHe);
            if (m->WS2CN != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2CN, m->residues[l].S2CNe);
            if (m->WS2CC != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2NHs*S2NHf, m->residues[l].S2CC, m->residues[l].S2CCe);
            fprintf(orderparams, "\n");
        } else if (m->model == MOD_SMF || m->model == MOD_SMFT) {
            Decimal S2 = (Decimal) m->residues[l].parameters[1];

            fprintf(orderparams, "%d\t", l+1);
            fprintf(orderparams, "%lf\t%lf\t%lf\t", S2, m->residues[l].S2NH, m->residues[l].S2NHe);
            if (m->WS2CH != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2, m->residues[l].S2CH, m->residues[l].S2CHe);
            if (m->WS2CN != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2, m->residues[l].S2CN, m->residues[l].S2CNe);
            if (m->WS2CC != 0) fprintf(orderparams, "%lf\t%lf\t%lf\t", S2, m->residues[l].S2CC, m->residues[l].S2CCe);
            fprintf(orderparams, "\n");
        }
    }
    fclose(gaf);
    fclose(orderparams);
    return 1;
}

int print_backcalcs(struct Model *m) {
    char filename[300];
    unsigned int l;
    for (l = m->proc_start; l < m->proc_end; l++) {
        sprintf(filename, "%s/backcalc_%d.dat", m->outputdir, l+1);
        back_calculate((m->residues[l].parameters), &(m->residues[l]), m, filename, m->params);
    }
    return 1;
}

int determine_residues(struct Model *m, int myid, int numprocs, unsigned int *start, unsigned int *end) {
    // we have m->n_residues split over numprocs.
    // so


    unsigned int res_per_proc = (m->n_residues / numprocs);
    if (m->n_residues % numprocs != 0)
        res_per_proc++;
    *start = res_per_proc * myid;
    *end = *start + res_per_proc;
    if (*end > m->n_residues)
        *end = m->n_residues;
    return 1;
}