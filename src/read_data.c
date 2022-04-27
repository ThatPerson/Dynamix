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
 * read_data.c
 *
 * Reads input data files and also provides utility to print out model file.
 */


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "datatypes.h"

int read_csa(struct Model *m, char *filename, int dt) {
    FILE *fp;
    char line[255];
    int len = 255;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        ERROR("%s not found.", filename);
        return -1;
    }
    int resid;
    Decimal s11, s22, s33;
    while (fgets(line, len, fp)) {
        if (line[0] == '%')
            continue; // comment
        int k = sscanf(line, "%d %lf %lf %lf\n", &resid, &s11, &s22, &s33);
        printf("%s -> %d, %lf, %lf, %lf\n", line, resid, s11, s22, s33);
        if (k != 4)
            ERROR("Reading %s failed.", line);
        resid = resid - 1; // 0 indexed in C

        if (s11 == -1 || s22 == -1 || s33 == -1)
            m->residues[resid].ignore = 1;
        switch (dt) {
            case DATA_CSISOC:
                m->residues[resid].csaC[0] = s11;
                m->residues[resid].csaC[1] = s22;
                m->residues[resid].csaC[2] = s33;
                break;
            case DATA_CSISON:
                m->residues[resid].csaN[0] = s11;
                m->residues[resid].csaN[1] = s22;
                m->residues[resid].csaN[2] = s33;
                break;
            default:
                break;
        }
    }
    fclose(fp);
    return 1;
}

/**
 * Reads residue specific data from datafile. \n
 * Data file being read should be set out in columns, with spaces between them.\n
 * \t Column one should contain the residue number as an integer, indexed from 1.\n
 * \t Column two should contain the data of note.\n
 * \t Column three should contain the error (or 0).\n
 * @note Comments may be inserted by preceeding a line with '%'\n
 * @param m
 *  Pointer to model
 * @param filename
 *  File to read
 * @param dt
 *  Data type (one of DATA_S2, DATA_CSISOC, DATA_CSISON)
 * @return 1 if successful, -1 else.
 */
int read_resid_data(struct Model *m, char *filename, int dt) {
    FILE *fp;
    char line[255];
    int len = 255;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        ERROR("%s not found.", filename);
        return -1;
    }
    int resid;
    Decimal val;
    Decimal err;

    while (fgets(line, len, fp)) {
        if (line[0] == '%')
            continue; // comment
        int k = sscanf(line, "%d %lf %lf\n", &resid, &val, &err);
        if (k != 3)
            ERROR("Reading %s failed.", line);
        resid = resid - 1; // 0 indexed in C
        switch (dt) {
            case DATA_S2NH:
                m->residues[resid].S2NH = val;
                m->residues[resid].S2NHe = err;
                if (val == -1) m->residues[resid].ignore = 1;
                if (err <= 0) m->residues[resid].S2NHe = 0.01;
                break;
            case DATA_S2CH:
                m->residues[resid].S2CH = val;
                m->residues[resid].S2CHe = err;
                if (val == -1 && (m->model == MOD_GAF || m->model == MOD_GAFT)) m->residues[resid].ignore = 1;
                if (err <= 0) m->residues[resid].S2CHe = 0.01;
                break;
            case DATA_S2CC:
                m->residues[resid].S2CC = val;
                m->residues[resid].S2CCe = err;
                if (val == -1 && (m->model == MOD_GAF || m->model == MOD_GAFT)) m->residues[resid].ignore = 1;
                if (err <= 0) m->residues[resid].S2CCe = 0.01;
                break;
            case DATA_S2CN:
                m->residues[resid].S2CN = val;
                m->residues[resid].S2CNe = err;
                if (val == -1 && (m->model == MOD_GAF || m->model == MOD_GAFT)) m->residues[resid].ignore = 1;
                if (err <= 0) m->residues[resid].S2CNe = 0.01;
                break;
            case DATA_CSISON:
                m->residues[resid].csisoN = val;
                /* 15N CSA parametrization according to isotropic chemical shift
                 * Wylie, B. J. et al.; Rienstra, C. M. Proc. Natl. Acad. Sci. U. S. A. 2011, 108, 16974.%lfrom Rienstra, PNAS 2011.:
                 */
                m->residues[resid].csaN[0] = 1.1283 * val + 93.77;
                m->residues[resid].csaN[1] = 1.0086 * val - 42.475;
                m->residues[resid].csaN[2] = 0.8631 * val - 51.295;
                break;
            case DATA_CSISOC:
                m->residues[resid].csisoC = val;
                /* 13C' CSA parametrization according to isotropic chemical shift from
                 * Wylie, B. J. et al.; Rienstra, C. M. J. Am. Chem. Soc. 2007, 129, 5318.*/
                m->residues[resid].csaC[0] = 0.24 * val + 200;
                m->residues[resid].csaC[1] = 2.82 * val - 305;
                m->residues[resid].csaC[2] = 96.5;
                //m->residues[resid].csaC[0] *= 1.08;
                //m->residues[resid].csaC[1] *= 1.08;
                //m->residues[resid].csaC[2] *= 1.08;

                break;
            default:
                break;
        }
    }
    fclose(fp);
    return 1;
}

/**
 * Reads peptide plane orientation data. Should be in a file with 3 spaced columns;\n
 * \t Column 1 should contain the residue number indexed from 1\n
 * \t Column 2 should contain the theta angle (in radians)\n
 * \t Column 3 should contain the phi angle (in radians)\n
 * @note Comments may be inserted by preceeding a line with '%'\n
 * @param m
 *  Pointer to model
 * @param filename
 *  Filename
 * @param orient
 *  One of OR_* denoting which orientations are being read.
 * @return 1 if successful, else -1.
 */
int read_pp(struct Model *m, char *filename, unsigned int orient) {
    FILE *fp;
    char line[255];
    int len = 255;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        ERROR("%s not found.", filename);
        return -1;
    }
    int resid;
    Decimal theta;
    Decimal phi;
    while (fgets(line, len, fp)) {
        if (line[0] == '%')
            continue; // comment
        int k = sscanf(line, "%d %lf %lf\n", &resid, &theta, &phi);
        if (k != 3)
            ERROR("Reading line '%s' failed.", line);
        resid = resid - 1; // index from 0
        m->residues[resid].orients[orient].theta = theta;
        m->residues[resid].orients[orient].phi = phi;
        calculate_Y2(&(m->residues[resid].orients[orient]));
    }
    fclose(fp);
    return 1;
}

/**
 * Reads relaxation data. Relaxation data should be set out with two components.\n
 * The file should begin with the keys 'FIELD', 'WR', 'W1', 'TYPE' denoting the form of the data, eg\n
 * \t FIELD = 700\n
 * \t WR = 50000\n
 * \t W1 = 0\n
 * \t TEMP = 300\n
 * \t TYPE = 15NR1\n
 * Each parameter will default to -1 so if one is not set this may cause errors (though a warning will be given). \n
 * TYPE may be one of [15NR1, 15NR1p, 13CR1, 13CR1p] currently, though it should be noted that 13C has not been tested.\n
 * @warning There must be a space on each size of the equals sign (eg, TEMP = 300 is good, TEMP=300 is bad)\n
 * Once this head section is complete, it should be followed by "#DATA".\n
 * After this, data should be split into three columns;\n
 * \t Column 1: residue number (starting from 1)\n
 * \t Column 2: Relaxation rate (in s-1)\n
 * \t Column 3: Relaxation error, two standard deviations (in s-1)\n
 * @note Comments may be inserted by preceeding a line with '%'\n
 * @param m
 *  Pointer to model
 * @param filename
 *  Filename
 * @return 1 if successful, else -1.
 */
int read_relaxation_data(struct Model *m, char *filename) {
    FILE *fp;
    char line[255];
    int len = 255;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        ERROR("%s not found.", filename);
        return -1;
    }
    int mode = 0;
    char key[255], val[255];
    unsigned int i;
    Decimal field = -1; // in MHz
    Decimal wr = -1, c_w1 = -1, c_wr; // in Hz
    Decimal w1 = -1; // in Hz
    int type = -1, compensate = NO_COMPENSATE;
    Decimal T = -1; // in Kelvin
    Decimal Gd = 0;
    int hydrogen = PROTONATED;
    int resid;
    int rel = -1;
    Decimal R, Re;
    while (fgets(line, len, fp)) {
        if (line[0] == '%' || strlen(line) < 2)
            continue; // comment
        if (strncmp(line, "#DATA", 5) == 0) {
            for (i = 0; i < m->n_residues; i++) {
                /* first check we have enough memory */
                if (m->residues[i].lim_relaxation - m->residues[i].n_relaxation < 5) {
                    printf("Oof\n");
                    m->residues[i].lim_relaxation += N_RELAXATION;
                    m->residues[i].relaxation = (struct Relaxation *) realloc(m->residues[i].relaxation,
                                                                              sizeof(struct Relaxation) *
                                                                              (m->residues[i].lim_relaxation));
                    if (m->residues[i].relaxation == NULL) {
                        ERROR("Pointer loss");
                        fclose(fp); // clean up on the way out
                        exit(-1);
                    }
                }
                rel = (int) m->residues[i].n_relaxation;
                m->residues[i].n_relaxation++;
                if (field == -1)
                    printf("Note: Field = 0 in %s\n", filename);
                if (wr == -1)
                    printf("Note: wr = 0 in %s\n", filename);
                if (w1 == -1)
                    printf("Note: w1 = 0 in %s\n", filename);
                if (type == -1)
                    printf("Note: type = 0 in %s\n", filename);
                if (T == -1)
                    printf("Note: T = 0 in %s\n", filename);

                m->residues[i].relaxation[rel].field = field;
                m->residues[i].relaxation[rel].wr = wr;
                m->residues[i].relaxation[rel].w1 = w1;
                m->residues[i].relaxation[rel].compensate = compensate;
                m->residues[i].relaxation[rel].compensate_w1 = ((compensate & COMPENSATE_W1) != 0 ? c_w1 : 0);
                m->residues[i].relaxation[rel].compensate_wr = ((compensate & COMPENSATE_WR) != 0 ? c_wr : 0);
                m->residues[i].relaxation[rel].type = type;
                m->residues[i].relaxation[rel].T = T;
                m->residues[i].relaxation[rel].Gd = Gd;
                m->residues[i].relaxation[rel].hydrogen = hydrogen; // protonation state.

            }
            mode = 1;
            continue;
        }
        if (mode == 0) {
            sscanf(line, "%255s = %255s\n", key, val);
            if (strcmp(key, "FIELD") == 0)
                field = atof(val);
            else if (strcmp(key, "WR") == 0)
                wr = atof(val);
            else if (strcmp(key, "W1") == 0)
                w1 = atof(val);
            else if (strcmp(key, "TEMP") == 0)
                T = atof(val);
            else if (strcmp(key, "GD") == 0)
                Gd = atof(val);
            else if (strcmp(key, "TYPE") == 0) {
                if (strcmp(val, "15NR1") == 0)
                    type = R_15NR1;
                else if (strcmp(val, "15NR1p") == 0)
                    type = R_15NR1p;
                else if (strcmp(val, "13CR1") == 0)
                    type = R_13CR1;
                else if (strcmp(val, "13CR1p") == 0)
                    type = R_13CR1p;
            } else if (strcmp(key, "PROTONATED") == 0) {
                hydrogen = PROTONATED;
            } else if (strcmp(key, "DEUTERATED") == 0) {
                hydrogen = DEUTERATED;
            } else if (strcmp(key, "COMPENSATEW1") == 0) {
                compensate += COMPENSATE_W1;
                c_w1 = atof(val);
            } else if (strcmp(key, "COMPENSATEWR") == 0) {
                compensate += COMPENSATE_WR;
                c_wr = atof(val);
            } else {
                printf("Parameter %s unknown.\n", key);
                fclose(fp);
                free_all(m);
                return -1;
            }
        } else if (mode == 1) {
            int k = sscanf(line, "%d %lf %lf\n", &resid, &R, &Re);
            if (k != 3) {
                ERROR("Error reading line '%s'", line);
            }/* else {
				printf("%s : %d, %lf %lf\n", line, resid, R, Re);
			}*/
            //printf("%d, %d, %lf, %lf\n", rel, resid, R, Re);
            // index from 0
            /*if (R == -1) {
                //printf("%s %d\n", line, m->residues[resid-1].n_relaxation);
                m->residues[resid-1].n_relaxation = m->residues[resid-1].n_relaxation - 1;
                if (m->residues[resid-1].n_relaxation < 0)
                    m->residues[resid-1].n_relaxation = 0;
                continue;
            } else {*/
            m->residues[resid - 1].relaxation[rel].R = R;
            m->residues[resid - 1].relaxation[rel].Rerror = Re; // note indexing from 0
            //}

        }
    }
    fclose(fp);
    return 1;
}

/**
 * Reads system file. System file should begin with key value pairs laid out as\n
 * \t KEY = VALUE\n
 * With spaces either side of the equals. Keys which may be used are as follows;\n
 * \t MODEL (one of SMF, SMFT, EMF, EMFT, GAF, GAFT)\n
 * \t S2DIP (file containing dipolar order parameters)\n
 * \t CSISON (file containing isotropic chemical shifts for N)\n
 * \t CSISOC (file containign isotropic chemical shifts for C)\n
 * \t N_RESIDUES (number of residues in protein)\n
 * \t OUTPUT (directory for output data)\n
 * \t N_ITER (number of iterations for main fitting) \n
 * \t IGNORE (may be repeated, once for each residue to ignore.)\n
 * \t N_ERROR_ITER (number of iterations for error calculation) \n
 * \t OR_* (file containing orientation data. May be one of OR_NH, OR_NC, OR_NCA, OR_NCSAxx, OR_NCSAyy, OR_NCSAzz, OR_CCAp, OR_CCAc, OR_CN, OR_CNH, OR_CCSAxx, OR_CCAyy, OR_CCAzz. Should be set out as in read_pp().)\n
 * \t NTHREADS (number of parallel threads to spawn)\n
 * This should then be followed by "#RELAXATION". After this, the files containing relaxation data (as set out in read_relaxation_data()) should be listed.
 * 
 * @note Comments may be inserted by preceeding a line with '%'
 * @warning There must be a space on each size of the equals sign (eg, TEMP = 300 is good, TEMP=300 is bad)\n
 * @param filename
 *  filename
 * @param m
 *  Pointer to model
 * @return 1 if successful, else -1.
 */
int read_system_file(char *filename, struct Model *m) {
    /* This function reads a *.dx file as laid out in the docs.
     * It reads in parameters to the struct m, and then calls
     * read_relaxation and relevant other functions to import
     * additional data.
     *
     * Function returns 1 if successful, -1 if failure.
     */

    FILE *fp;
    char line[255];
    int len = 255;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        ERROR("%s not found.", filename);
        return -1;
    }
    int mode = 0;
    char key[255];
    char val[255];
    char s2nh[255] = "", s2ch[255] = "", s2cc[255] = "", s2cn[255] = "";
    char csisoN[255] = "";
    char csisoC[255] = "";
    char csaN[255] = "";
    char csaC[255] = "";
    int n_resid = -1;
    char pp_orient[N_OR][255];
    unsigned int i;
    for (i = 0; i < N_OR; i++) {
        strcpy(pp_orient[i], "");
    }
    strcpy(m->outputdir, "./");
    unsigned int to_ignore[128];
    int ig = 0;

    /* Initialise */
    m->n_error_iter = 0;
    m->n_bc_iter = 1;
    m->model = MOD_UNDEFINED;
    m->or_variation = INVARIANT_A;
    m->n_residues = 0;
    m->params = 0;
    m->nthreads = 1;
    m->cn_ratio = CNRATIO_OFF;
    m->global = LOCAL;
    m->gd_mod = GD_NO_MOD;

    m->n_anneal_iter = 1;
    m->n_nm_iter = 1;

    m->WS2NH = 1;
    m->WS2CH = 1;
    m->WS2CN = 1;
    m->WS2CC = 1;

    m->t_error = 0;

    m->anneal_restart = 0.01;
    m->anneal_wobb = 0.05;
    m->anneal_temp = 6000;
    m->anneal_therm = 1.1;
    m->ultrafast = DISABLED;
    m->microsecond = DISABLED;
    m->GDS2 = -1;
    m->GDtaur = -1;
    m->OValpha = -1;
    m->OVbeta = -1;
    m->OVgamma = -1;
    m->UFS2 = -1;
    while (fgets(line, len, fp)) {
        if (line[0] == '%')
            continue; // comment
        if (strncmp(line, "#RELAXATION", 11) == 0) {
            if (n_resid == -1) {
                ERROR("Number of residues must be defined.");
                return -1;
            }

            if (m->gd_mod == GD_MOD && m->ultrafast == ENABLED && m->or_variation == VARIANT_A) {
                m->GDS2 = m->params - 7;
                m->GDtaur = m->params - 6;
                m->OValpha = m->params - 5;
                m->OVbeta = m->params - 4;
                m->OVgamma = m->params - 3;
                m->UFtau_uf = m->params - 2;
                m->UFS2 = m->params - 1;
            } else if (m->gd_mod == GD_MOD && m->ultrafast == ENABLED) {
                m->GDS2 = m->params - 4;
                m->GDtaur = m->params - 3;
                m->UFtau_uf = m->params - 2;
                m->UFS2 = m->params - 1;
            } else if (m->gd_mod == GD_MOD && m->or_variation == VARIANT_A) {
                m->GDS2 = m->params - 5;
                m->GDtaur = m->params - 4;
                m->OValpha = m->params - 3;
                m->OVbeta = m->params - 2;
                m->OVgamma = m->params - 1;
            } else if (m->ultrafast == ENABLED && m->or_variation == VARIANT_A) {
                m->OValpha = m->params - 5;
                m->OVbeta = m->params - 4;
                m->OVgamma = m->params - 3;
                m->UFtau_uf = m->params - 2;
                m->UFS2 = m->params - 1;
            } else if (m->gd_mod == GD_MOD) {
                m->GDS2 = m->params - 2;
                m->GDtaur = m->params - 1;
            } else if (m->ultrafast == ENABLED) {
                m->UFtau_uf = m->params - 2;
                m->UFS2 = m->params - 1;
            } else if (m->or_variation == VARIANT_A) {
                m->OValpha = m->params - 3;
                m->OVbeta = m->params - 2;
                m->OVgamma = m->params - 1;
            }

            m->residues = (struct Residue *) malloc(sizeof(struct Residue) * (long unsigned int) n_resid);
            for (i = 0; i < (unsigned int) n_resid; i++) {
                m->residues[i].relaxation = (struct Relaxation *) malloc(sizeof(struct Relaxation) * N_RELAXATION);
                m->residues[i].lim_relaxation = N_RELAXATION;
                m->residues[i].n_relaxation = 0;
                m->residues[i].ignore = 0;
            }

            for (i = 0; i < (unsigned int) ig; i++) {
                m->residues[to_ignore[i] - 1].ignore = 1;
            }

            int t = 0, q = 0; // t is general counter, q is for GAF models only
            if (strcmp(s2nh, "") != 0)
                t += read_resid_data(m, s2nh, DATA_S2NH);
            else {
                m->WS2NH = 0;
            }

            if (strcmp(s2cc, "") != 0)
                q += read_resid_data(m, s2cc, DATA_S2CC);
            else {
                q++;
                m->WS2CC = 0;
            }

            if (strcmp(s2ch, "") != 0)
                q += read_resid_data(m, s2ch, DATA_S2CH);
            else {
                q++;
                m->WS2CH = 0;
            }

            if (strcmp(s2cn, "") != 0)
                q += read_resid_data(m, s2cn, DATA_S2CN);
            else {
                q++;
                m->WS2CN = 0;
            }


            if (strcmp(csisoN, "") != 0)
                t += read_resid_data(m, csisoN, DATA_CSISON);
            if (strcmp(csaN, "") != 0)
                read_csa(m, csaN, DATA_CSISON);


            if (strcmp(csisoC, "") != 0)
                t += read_resid_data(m, csisoC, DATA_CSISOC);
            if (strcmp(csaC, "") != 0)
                read_csa(m, csaC, DATA_CSISOC);

            if (t != 3) {
                ERROR("error reading one of S2*, csisoN, csisoC");
                return -1;
            }
            if (q != 3 && (m->model == MOD_GAF || m->model == MOD_GAFT)) {
                ERROR("error reading one of S2CH, S2CN, S2CC");
                ERROR("The weight of those not present was set to 0.");
            }

            t = 0;
            for (i = 0; i < N_OR; i++) {
                if (strcmp(pp_orient[i], "") != 0)
                    t += read_pp(m, pp_orient[i], i);
                else
                    t++;
            }
            if (t != N_OR) {
                ERROR("error in reading orientations");
                return -1;
            }
            mode = 1;
            continue;
        }
        if (mode == 0) {
            if (line[0] == '#')
                continue;
            sscanf(line, "%255s = %255s\n", key, val);
            if (strcmp(key, "MODEL") == 0) {
                if (strcmp(val, "SMF") == 0) {
                    m->params = 2;
                    m->model = MOD_SMF;
                } else if (strcmp(val, "EMF") == 0) {
                    m->params = 3;
                    m->model = MOD_EMF;
                } else if (strcmp(val, "EMFT") == 0) {
                    m->params = 5;
                    m->model = MOD_EMFT;
                } else if (strcmp(val, "SMFT") == 0) {
                    m->params = 3;
                    m->model = MOD_SMFT;
                } else if (strcmp(val, "GAF") == 0) {
                    m->params = 8;
                    m->model = MOD_GAF;
                } else if (strcmp(val, "GAFT") == 0) {
                    m->params = 10;
                    m->model = MOD_GAFT;
                } else if (strcmp(val, "DEMF") == 0) {
                    m->params = 4;
                    m->model = MOD_DEMF;
                } else if (strcmp(val, "DEMFT") == 0) {
                    m->params = 6;
                    m->model = MOD_DEMFT;
                } else if (strcmp(val, "RDEMFT") == 0) {
                    m->params = 8;
                    m->model = MOD_RDEMFT;
                } else if (strcmp(val, "SDEMFT") == 0) {
                    m->params = 8;
                    m->model = MOD_SDEMFT;
                } else if (strcmp(val, "SDEMF") == 0) {
                    m->params = 6;
                    m->model = MOD_SDEMF;
                } else if (strcmp(val, "EGAF") == 0) {
                    m->params = 6;
                    m->model = MOD_EGAF;
                } else if (strcmp(val, "EGAFT") == 0) {
                    m->params = 8;
                    m->model = MOD_EGAFT;
                } else if (strcmp(val, "BGAF") == 0) {
                    m->params = 4;
                    m->model = MOD_BGAF;
                } else if (strcmp(val, "BGAFT") == 0) {
                    m->params = 5;
                    m->model = MOD_BGAFT;
                } else if (strcmp(val, "BAIMF") == 0) {
                    m->params = 4;
                    m->model = MOD_BAIMF;
                } else if (strcmp(val, "BAIMFT") == 0) {
                    m->params = 5;
                    m->model = MOD_BAIMFT;
                } else if (strcmp(val, "AIMF") == 0) {
                    m->params = 8;
                    m->model = MOD_AIMF;
                } else if (strcmp(val, "AIMFT") == 0) {
                    m->params = 10;
                    m->model = MOD_AIMFT;
                } else if (strcmp(val, "EAIMF") == 0) {
                    m->params = 6;
                    m->model = MOD_EAIMF;
                } else if (strcmp(val, "EAIMFT") == 0) {
                    m->params = 8;
                    m->model = MOD_EAIMFT;
                } else if (strcmp(val, "BGF") == 0) {
                    m->params = 4;
                    m->model = MOD_BGF;
                } else if (strcmp(val, "BGFT") == 0) {
                    m->params = 6;
                    m->model = MOD_BGFT;
                } else {
                    printf("Model %s unknown.\n", val);
                    return -1;
                }
            } else if (strcmp(key, "S2NH") == 0) {
                strcpy(s2nh, val);
            } else if (strcmp(key, "S2CH") == 0) {
                strcpy(s2ch, val);
            } else if (strcmp(key, "S2CC") == 0) {
                strcpy(s2cc, val);
            } else if (strcmp(key, "S2CN") == 0) {
                strcpy(s2cn, val);
            } else if (strcmp(key, "W_S2NH") == 0) {
                m->WS2NH = atoi(val);
            } else if (strcmp(key, "W_S2CH") == 0) {
                m->WS2CH = atoi(val);
            } else if (strcmp(key, "W_S2CN") == 0) {
                m->WS2CN = atoi(val);
            } else if (strcmp(key, "W_S2CC") == 0) {
                m->WS2CC = atoi(val);
            } else if (strcmp(key, "OR_VARY") == 0) {
                if (m->or_variation == VARIANT_A)
                    continue;
                m->params += 3;
                m->or_variation = VARIANT_A;
            } else if (strcmp(key, "GDMOD") == 0) {
                if (m->gd_mod == GD_MOD)
                    continue;
                m->params += 2;
                m->gd_mod = GD_MOD;
            } else if (strcmp(key, "ULTRAFAST") == 0) {
                if (m->ultrafast == ENABLED)
                    continue;
                m->params += 2;
                m->ultrafast = ENABLED;
             //   m->microsecond = ENABLED;
            } else if (strcmp(key, "MICROSECOND") == 0) {
                m->microsecond = (atoi(val) == 1 ? ENABLED : DISABLED);
            } else if (strcmp(key, "CSISON") == 0) {
                strcpy(csisoN, val);
            } else if (strcmp(key, "CSISOC") == 0) {
                strcpy(csisoC, val);
            } else if (strcmp(key, "CSAN") == 0) {
                strcpy(csaN, val);
            } else if (strcmp(key, "CSAC") == 0) {
                strcpy(csaC, val);
            } else if (strcmp(key, "N_RESIDUES") == 0) {
                m->n_residues = (unsigned int) atoi(val);
                n_resid = atoi(val);
            } else if (strcmp(key, "N_NM_ITER") == 0) {
                m->n_nm_iter = (unsigned int) atoi(val);
            } else if (strcmp(key, "N_ANNEAL_ITER") == 0) {
                m->n_anneal_iter = (unsigned int) atoi(val);
            } else if (strcmp(key, "N_BC_ITER") == 0) {
                m->n_bc_iter = (unsigned int) atoi(val);
            } else if (strcmp(key, "ANNEAL_TEMP") == 0) {
                m->anneal_temp = atof(val);
            } else if (strcmp(key, "ANNEAL_WOBB") == 0) {
                m->anneal_wobb = atof(val);
            } else if (strcmp(key, "ANNEAL_THERM") == 0) {
                m->anneal_therm = atof(val);
            } else if (strcmp(key, "ANNEAL_RESTART") == 0) {
                m->anneal_restart = atof(val);
            } else if (strcmp(key, "OUTPUT") == 0) {
                strcpy(m->outputdir, val);
            } else if (strcmp(key, "IGNORE") == 0) {
                to_ignore[ig] = (unsigned int) atoi(val);
                ig++;
            } else if (strcmp(key, "CNCOMP") == 0) {
                m->cn_ratio = CNRATIO_ON;
            } else if (strcmp(key, "TERROR") == 0) {
                m->t_error = atof(val);
            } else if (strcmp(key, "GLOBAL") == 0) {
                m->global = GLOBAL;
            } else if (strcmp(key, "N_ERROR_ITER") == 0) {
                m->n_error_iter = (unsigned int) atoi(val);
            } else if (strcmp(key, "OR_NH") == 0) {
                strcpy(pp_orient[OR_NH], val);
            } else if (strcmp(key, "OR_NC") == 0) {
                strcpy(pp_orient[OR_NC], val);
            } else if (strcmp(key, "OR_NCA") == 0) {
                strcpy(pp_orient[OR_NCA], val);
            } else if (strcmp(key, "OR_NCSAxx") == 0) {
                strcpy(pp_orient[OR_NCSAxx], val);
            } else if (strcmp(key, "OR_NCSAyy") == 0) {
                strcpy(pp_orient[OR_NCSAyy], val);
            } else if (strcmp(key, "OR_NCSAzz") == 0) {
                strcpy(pp_orient[OR_NCSAzz], val);
            } else if (strcmp(key, "OR_CCAp") == 0) {
                strcpy(pp_orient[OR_CCAp], val);
            } else if (strcmp(key, "OR_CCAc") == 0) {
                strcpy(pp_orient[OR_CCAc], val);
            } else if (strcmp(key, "OR_CN") == 0) {
                strcpy(pp_orient[OR_CN], val);
            } else if (strcmp(key, "OR_CNH") == 0) {
                strcpy(pp_orient[OR_CNH], val);
            } else if (strcmp(key, "OR_CH") == 0) {
                strcpy(pp_orient[OR_CH], val);
            } else if (strcmp(key, "OR_CCSAxx") == 0) {
                strcpy(pp_orient[OR_CCSAxx], val);
            } else if (strcmp(key, "OR_CCSAyy") == 0) {
                strcpy(pp_orient[OR_CCSAyy], val);
            } else if (strcmp(key, "OR_CCSAzz") == 0) {
                strcpy(pp_orient[OR_CCSAzz], val);
            } else if (strcmp(key, "NTHREADS") == 0) {
                m->nthreads = (unsigned int) atoi(val);
                if (m->nthreads <= 0) {
                    printf("Number of threads must be greater than 0.\n");
                    exit(-1);
                }
            } else {
                ERROR("Unrecognized argument %s (%s)\n", key, val);
            }

            //printf("%s: %s\n", key, val);
        } else if (mode == 1) {
            strcpy(key, "IGNORE");
            sscanf(line, "%255s\n", key);
            //printf("Read file: %s\n", key);
            //int read_relaxation_data(struct Model *m, char *filename) {
            if (strcmp(key, "IGNORE") != 0)
                read_relaxation_data(m, key);
        }
    }
    unsigned int k;
    int c_count, n_count;
    for (i = 0; i < m->n_residues; i++) {
        //if (m->residues[i].n_relaxation < 5)
        //	m->residues[i].ignore = 1;
        if (m->cn_ratio == CNRATIO_OFF)
            continue;
        c_count = 0;
        n_count = 0;
        for (k = 0; k < m->residues[i].n_relaxation; k++) {
            if (m->residues[i].relaxation[k].type == R_13CR1 || m->residues[i].relaxation[k].type == R_13CR1p)
                c_count++;
            else if (m->residues[i].relaxation[k].type == R_15NR1 || m->residues[i].relaxation[k].type == R_15NR1p)
                n_count++;
        }
        if (n_count == 0 || c_count == 0) {
            printf("Residue %d: Only one nuclei measured. CN ratio set to 1.\n", i + 1);
            m->residues[i].cn = 1;
        } else {
            m->residues[i].cn = ((Decimal) n_count) / ((Decimal) c_count);
            //printf("%d, %d -> %lf\n", c_count, n_count, m->residues[i].cn);
        }
    }
    fclose(fp);

    /*Check requirements */
    if (m->n_nm_iter == 0 || m->n_anneal_iter == 0) {
        ERROR("Please provide N_NM_ITER and N_ANNEAL_ITER");
        return -1;
    }
    if (m->model == MOD_UNDEFINED) {
        ERROR("Please provide MODEL");
        return -1;
    }
    if (m->n_residues == 0) {
        ERROR("Please provide residue count");
        return -1;
    }

    return 1;


    /*	int max_func_evals, max_iter;
    int model;
    struct Residue * residues;
    int n_residues;*/

}

int print_system(struct Model *m, char *filename) {
    FILE *fp;
    fp = fopen(filename, "w");
    if (fp == NULL) {
        ERROR("%s not found.", filename);
        return -1;
    }


    //fprintf(fp, "MFE: %d\nMI:  %d\n", m->max_func_evals, m->max_iter);
    fprintf(fp, "Model: %d\nN_Residues: %d\n", m->model, m->n_residues);
    fprintf(fp, "Params: %d\nN threads: %d\n", m->params, m->nthreads);
    fprintf(fp, "Error Mode: %s\n", (m->error_mode == 1) ? "ON" : "OFF");
    fprintf(fp, "Orientation Variation: %s\n", (m->or_variation == VARIANT_A) ? "ON" : "OFF");
    if (m->or_variation == VARIANT_A)
        fprintf(fp, "\tOV alpha: %d\n\tOV beta: %d\n\tOV gamma: %d\n", m->OValpha, m->OVbeta, m->OVgamma);
    fprintf(fp, "Ultrafast: %s\n", (m->ultrafast == ENABLED) ? "ON" : "OFF");
    if (m->ultrafast == ENABLED)
        fprintf(fp, "\tUF S2: %d\n\tUF tuf: %d\n", m->UFS2, m->UFtau_uf);
    fprintf(fp, "us Timescale: %s\n", (m->microsecond == ENABLED) ? "ON" : "OFF");
    if (m->gd_mod == GD_MOD)
        fprintf(fp, "\tGD S2: %d\n\tGD tr: %d\n", m->GDS2, m->GDtaur);
    fprintf(fp, "Gd PRE: %s\n", (m->gd_mod == GD_MOD)?"ON":"OFF");
    fprintf(fp, "C/N Ratio Compensation: %s\n", (m->cn_ratio == CNRATIO_ON) ? "ON" : "OFF");
    fprintf(fp, "Global: %s\n", (m->global == GLOBAL) ? "ON" : "OFF");

    unsigned int i, j;
    struct Relaxation *r;
    for (i = 0; i < m->n_residues; i++) {
        fprintf(fp, "=== Residue %d ===\n", i + 1);
        if (m->residues[i].ignore == 1)
            fprintf(fp, "--- IGNORING ---\n");

        fprintf(fp, "\tS2NH = %lf (%d) \n\tS2CH = %lf (%d) \n\tS2CC = %lf (%d) \n\tS2CN = %lf (%d) \n",
                m->residues[i].S2NH, m->WS2NH, m->residues[i].S2CH, m->WS2CH, m->residues[i].S2CC, m->WS2CC,
                m->residues[i].S2CN, m->WS2CN);
        fprintf(fp, "\tCSISON = %lf\n\tCSISOC = %lf\n", m->residues[i].csisoN, m->residues[i].csisoC);
        fprintf(fp, "\tCSAN: [%lf, %lf, %lf]\n", m->residues[i].csaN[0], m->residues[i].csaN[1],
                m->residues[i].csaN[2]);
        Decimal Qcsiso = m->residues[i].csisoN;
        Decimal Qred_aniso = m->residues[i].csaN[0] - Qcsiso;
        Decimal Qeta = (m->residues[i].csaN[1] - m->residues[i].csaN[2]) / Qred_aniso;
        fprintf(fp, "\t\tCSA: d %lf, n %lf\n", Qred_aniso, Qeta);
        // z y x
        fprintf(fp, "\tCSAC: [%lf, %lf, %lf]\n", m->residues[i].csaC[0], m->residues[i].csaC[1],
                m->residues[i].csaC[2]);
        Decimal Ccsiso = m->residues[i].csisoC;
        Decimal Cred_aniso = m->residues[i].csaC[0] - Ccsiso;
        Decimal Ceta = (m->residues[i].csaC[1] - m->residues[i].csaC[2]) / Cred_aniso;
        fprintf(fp, "\t\tCSA: d %lf, n %lf\n", Cred_aniso, Ceta);
        if (m->cn_ratio == CNRATIO_ON)
            fprintf(fp, "\tCNRATIO: %lf\n", m->residues[i].cn);
        else
            fprintf(fp, "\tCNRATIO: OFF\n");
        fprintf(fp, "\tOrients;\n");
        for (j = 0; j < N_OR; j++) {
            //if (m->residues[i].orients[j] != NULL)
            fprintf(fp, "\t\t%d: %lf, %lf\n", j, m->residues[i].orients[j].theta, m->residues[i].orients[j].phi);
        }
        fprintf(fp, "\tRelaxation Constraints: %d\n", m->residues[i].n_relaxation);
        for (j = 0; j < m->residues[i].n_relaxation; j++) {
            r = &(m->residues[i].relaxation[j]);
            fprintf(fp, "\t\t%d: %d [%lf, %lf, %lf, %lf] %lf +- %lf %s\n", j, r->type, r->field, r->wr, r->w1, r->T, r->R,
                    r->Rerror, (r->hydrogen == PROTONATED ? "PROTONATED":"DEUTERATED"));
        }
    }
    fclose(fp);
    return 1;

}
