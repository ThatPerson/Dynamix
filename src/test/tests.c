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
 * tests.c
 *
 * Unit testing using CMocka library.
 */

#include "tests.h"
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <complex.h>

#include <stddef.h>
#include <setjmp.h>
#include <stdint.h>
#include <cmocka.h> // CMocka unit testing

#include "../datatypes.h"
#include "../read_data.h"
#include "../models/smf.h"
#include "../models/emf.h"
#include "../models/gaf.h"
#include "../models/egaf.h"
#include "../models/aimf.h"
#include "../chisq.h"
#include "../crosen.h" // implementation of Nelder-Mead simplex algorithm
#include "anneal.h" // simulated annealing algorithm

#include "../errors.h"
#include "../global_gaf.h"
#include "../verification.h"
#include "../runners.h"

static void test_determine_residues(void **state) {
    (void) state;
    struct Model m;
    unsigned int start, end;
    m.n_residues = 56;
    int k = determine_residues(&m, 2, 4, &start, &end);
    assert_int_equal(k, 1);
    assert_int_equal(start, 28);
    assert_int_equal(end, 42);
    k = determine_residues(&m, 5, 9, &start, &end);
    assert_int_equal(k, 1);
    assert_int_equal(start, 35);
    assert_int_equal(end, 42);
    k = determine_residues(&m, 3, 4, &start, &end);
    assert_int_equal(k, 1);
    assert_int_equal(start, 42);
    assert_int_equal(end,56 );
    k = determine_residues(&m, 2, 3, &start, &end);
    assert_int_equal(k, 1);
    assert_int_equal(start, 38);
    assert_int_equal(end, 56);
}

static void test_temp_tau(void ** state) {
    (void) state;
    Decimal tau0 = pow(10, -4);
    Decimal epsilon = 0.0001;
    assert_float_equal(temp_tau(tau0, 0, 200), tau0, epsilon);
    assert_float_equal(temp_tau(tau0, 10, 200), 0.0001, epsilon); // small activation energy
    assert_float_equal(temp_tau(tau0, 10000, 200), 0.040896, epsilon);

}

static void test_dwig(void **state) {
    (void) state;
    Decimal beta = M_PI_2;
    Decimal epsilon = 0.0001;
    Decimal test_dwig[5][5];
    initialise_dwig(beta, test_dwig);
    assert_float_equal(test_dwig[2][0], 0.612372, epsilon);
    assert_float_equal(test_dwig[0][2], 0.612372, epsilon);
    assert_float_equal(test_dwig[2][2], -0.5, epsilon);
    assert_float_equal(test_dwig[4][4], 0.25, epsilon);
    assert_float_equal(test_dwig[1][4], 0.5, epsilon);
    beta = 0.64;
    initialise_dwig(beta, test_dwig);
    assert_float_equal(test_dwig[2][0], 0.218392, epsilon);
    assert_float_equal(test_dwig[0][2], 0.218398, epsilon);
    assert_float_equal(test_dwig[2][2], 0.465036, epsilon);
    assert_float_equal(test_dwig[4][4], 0.811887, epsilon);
    assert_float_equal(test_dwig[1][4], 0.0590938, epsilon);

}

static void test_sphericals(void **state) {
    (void) state;
    struct Orient test;
    test.theta = 3.82;
    test.phi = 2.42;
    calculate_Y2(&test);
    Decimal epsilon = 0.0001;
    assert_float_equal(creal(test.Y2[0]), 0.0193601, epsilon);
    assert_float_equal(cimag(test.Y2[0]), 0.150887, epsilon);
    assert_float_equal(creal(test.Y2[2]), 0.258157, epsilon);
    assert_float_equal(creal(test.Y2[3]), 0.283383, epsilon);
    assert_float_equal(cimag(test.Y2[4]), -0.150887, epsilon);
    unsigned int i;
    for (i = 0; i < 5; i++) {
        assert_float_equal(creal(test.Y2[i]), creal(test.Y2c[i]), epsilon);
        assert_float_equal(cimag(test.Y2[i]), -cimag(test.Y2c[i]), epsilon);
    }
}

static void test_rotations(void **state) {
    (void) state;
    struct Orient CACA, CACAp; // ca ca axis
    struct Orient pp, ppp; // plane perpendicula
    struct Orient NHv, NHvp; // roughly parallel to N-H
    unsigned int i;
    Decimal epsilon = 0.0001;
    CACA.theta = 0;
    CACA.phi = 0;
    pp.theta = HALF_PI;
    pp.phi = HALF_PI;
    NHv.theta = HALF_PI ;
    NHv.phi = 0;
    CACAp.theta = 0;
    CACAp.phi = 0;
    ppp.theta = HALF_PI;
    ppp.phi = HALF_PI;
    NHvp.theta = HALF_PI ;
    NHvp.phi = 0;

    calculate_Y2(&CACA);
    calculate_Y2(&pp);
    calculate_Y2(&NHv);
    calculate_Y2(&CACAp);
    calculate_Y2(&ppp);
    calculate_Y2(&NHvp);

    rotate_Y2(&CACA, HALF_PI, 0, 0);
    rotate_Y2(&pp, HALF_PI, 0, 0);
    rotate_Y2(&NHv, HALF_PI, 0, 0);

    for (i = 0; i < 5; i++) {
        assert_float_equal(creal(CACA.Y2[i]), creal(CACAp.Y2[i]), epsilon);
        assert_float_equal(cimag(CACA.Y2[i]), cimag(CACAp.Y2[i]), epsilon);
    }

    calculate_Y2(&CACA);
    calculate_Y2(&pp);
    calculate_Y2(&NHv);

    rotate_Y2(&CACA, 0, HALF_PI, 0);
    rotate_Y2(&pp, 0, HALF_PI, 0);
    rotate_Y2(&NHv, 0, HALF_PI, 0);

    for (i = 0; i < 5; i++) {
        assert_float_equal(creal(pp.Y2[i]), creal(ppp.Y2[i]), epsilon);
        assert_float_equal(cimag(pp.Y2[i]), cimag(ppp.Y2[i]), epsilon);
    }

}

static void single_gaf(Decimal sa, Decimal sb, Decimal sg, Decimal Atheta, Decimal Aphi, Decimal Btheta, Decimal Bphi, Decimal SAAr, Decimal SBBr, Decimal SABr) {
    Decimal epsilon = 0.0001;
    struct Orient A, B;
    A.theta = Atheta;
    A.phi = Aphi;
    B.theta = Btheta;
    B.phi = Bphi;
    calculate_Y2(&A);
    calculate_Y2(&B);
    struct Orient *set1[] = {&A, &B, &A};
    struct Orient *set2[] = {&A, &B, &B};
    Decimal SAA, SBB, SAB;
    Decimal *S2[] = {&SAA, &SBB, &SAB};
    Decimal sigs[] = {sa, sb, sg};
    GAF_S2(sigs, set1, set2, S2, 3, MODE_REAL);
    int i;
    assert_float_equal(SAA, SAAr, epsilon);
    assert_float_equal(SBB, SBBr, epsilon);
    assert_float_equal(SAB, SABr, epsilon);

    return;
}

static void test_gaf(void **state) {
    (void) state;
    initialise_dwig(HALF_PI, Dwig);
    Decimal epsilon = 0.0001;
    Decimal sigs[3] = {0.08, 0.41, 0.80};

    struct Orient A, B;
    A.theta = 2.03;
    A.phi = 1.64;
    B.theta = 4.07;
    B.phi = 1.84;
    calculate_Y2(&A);
    calculate_Y2(&B);

    Decimal SAA = 1., SBB = -1., SAB;
    struct Orient *set1[] = {&A, &B, &A};
    struct Orient *set2[] = {&A, &B, &B};
    Decimal *S2[] = {&SAA, &SBB, &SAB};

    GAF_S2(sigs, set1, set2, S2, 3, MODE_REAL);

    assert_float_equal(SAA, .279796, epsilon);
    assert_float_equal(SBB, .311143, epsilon);
    assert_float_equal(SAB, -0.218714, epsilon);

  //  static void single_gaf(Decimal sa, Decimal sb, Decimal sg, Decimal Atheta, Decimal Aphi, Decimal Btheta, Decimal Bphi, Decimal SAAr, Decimal SBBr, Decimal SABr)
    single_gaf(0.08, 0.41, 0.80, 2.03, 1.64, 4.07, 1.84, 0.279796, 0.311143,-0.218714);
    single_gaf(0.88, 0.02, 0.08, 4.79, M_PI, 3.35, M_PI, 0.970387, 0.277148, -0.448971);
    single_gaf(0.12, 0.98, 0.60, 5.62, HALF_PI, 1.34, HALF_PI, 0.187247, 0.332697, -0.121306);
}

static void test_statistics(void **state) {
    (void) state;
    int i;
    Decimal mean = (rand()%1000)/1000.;
    Decimal std = (rand()%1000)/1000.;
    Decimal calc_mean, calc_std;

    Decimal rnds[200000];
    for (i = 0; i < 200000; i++) {
        rnds[i] = (Decimal) norm_rand(mean, std);
    }
    calc_statistics(rnds, 200000, &calc_mean, &calc_std);
    assert_float_equal(mean, calc_mean, 0.01);
    assert_float_equal(std, calc_std, 0.01);

}

static void test_crosen_backcalc(void **state) {
    (void) state;
    struct Model m;
    //Decimal min = simplex(optimize_chisq, opts, 1.0e-16, 1, resid, m);

    // initialize a basic relaxation set from known parameters. Then backcalculate rates, and then fit those. Assert that the input is equal to the output

    m.model = MOD_DEMFT;
    m.params = 6;
    m.n_residues = 1;
    m.rdc = RDC_OFF;
    m.nthreads = 1;
    m.global = LOCAL;
    m.cn_ratio = CNRATIO_OFF;
    m.or_variation = INVARIANT_A;
    m.error_mode = 0;
    m.proc_start = 0;
    m.proc_end = 1;
    m.myid = 0;
    m.numprocs = 1;
    m.WS2NH = 100;
    m.WS2CC = 0;
    m.WS2CN = 0;
    m.WS2CH = 0;

    /* Taken from GB1 Residue 44 */

    m.residues = (struct Residue *) malloc(sizeof(struct Residue) * 1);
    if (m.residues == NULL) goto fail;
    m.residues[0].S2NH = 0.8;
    m.residues[0].S2CH = 0.832;
    m.residues[0].S2CC = 0.975;
    m.residues[0].S2CN = 0.906;
    m.residues[0].csisoN = 109.21;
    m.residues[0].csisoC = 176.72;
    m.residues[0].cn = 1.047619;
    m.residues[0].S2NHe = 0.01;
    m.residues[0].S2CNe = 0.05;
    m.residues[0].S2CHe = 0.05;
    m.residues[0].S2CCe = 0.05;
    m.residues[0].ignore = 0;
    m.residues[0].parameters = (Decimal *) malloc(sizeof(Decimal) * 6);
    if (m.residues[0].parameters == NULL) goto fail;
    unsigned int i;
    Decimal parms[] = {1e-04, 0.84211, 1e-03, 0.95, 3e+04, 4e+03};
    for (i = 0 ; i < 6; i++)
        m.residues[0].parameters[i] = parms[i];

    int ignore = 0;

    Decimal fields[] = {600, 800, 1000};
    Decimal spin_rates[] = {50000, 75000, 100000};
    Decimal nut_freq[] = {10000, 20000};
    Decimal temps[] = {250, 300, 350};
    unsigned int n_f, n_sr, n_nf, n_T, n_t;
    unsigned int N_f = 3, N_sr = 3, N_nf = 2, N_T = 3, N_t = 4;
    unsigned int N_rates = N_f * N_sr * N_nf * N_T * N_t;
    m.residues[0].relaxation = (struct Relaxation *) malloc(sizeof(struct Relaxation) * N_rates);
    m.residues[0].lim_relaxation = N_rates;
    m.residues[0].n_relaxation = N_rates;
    if (m.residues[0].relaxation == NULL)
        goto fail;
    unsigned int k = 0;
    struct Relaxation *r;
    struct Residue *resid = &(m.residues[0]);
    for (n_f = 0; n_f < N_f; n_f++) {
        for (n_sr = 0; n_sr < N_sr; n_sr++) {
            for (n_nf = 0; n_nf < N_nf; n_nf++) {
                for (n_T = 0; n_T < N_T; n_T++) {
                    for (n_t = 0; n_t < N_t; n_t++) {
                        r = &(m.residues[0].relaxation[k]);
                        r->field = fields[n_f];
                        r->wr = spin_rates[n_sr];
                        r->w1 = nut_freq[n_nf];
                        r->T = temps[n_T];
                        r->type = n_t;
                        r->R = 1;
                        r->R = back_calc(resid->parameters, resid, r, &m, &ignore);
                        r->Rerror = 0.2 * r->R;
                        k++;
                    }
                }
            }
        }
    }

    assert_int_equal(k, N_rates);
    Decimal opts[6]; // = {5e-04, 0.82, 0.43e-03, 0.99, 2e+03, 5e+04};
    Decimal min = 1000000;



    const Decimal minv[] = {0.00001,
                            resid->S2NH,
                            0.0001,
                            resid->S2NH,
                            0,
                            0};
    const Decimal maxv[] = {0.001,
                            1,
                            0.01,
                            1,
                            60000,
                            60000};



    min = anneal(optimize_chisq, opts, minv, maxv, 6, 6000, 10000, 0.05, 1.1, 0.01, NULL, MODE_RANDOM_RESTART+MODE_RANDOM_START, resid, &m);
    min = simplex(optimize_chisq, opts, 1.0e-16, 1, resid, &m);


    assert_float_equal(min, 0, 3);

    Decimal temp_R;
    for (k = 0; k < N_rates; k++) {
        r = &(m.residues[0].relaxation[k]);
        temp_R = back_calc(opts, resid, r, &m, &ignore);
        assert_float_equal(temp_R, r->R, 2*r->Rerror);
    }

    free(m.residues[0].relaxation);
    free(m.residues[0].parameters);
    free(m.residues);
    return;
    fail:
        free(m.residues[0].relaxation);
        free(m.residues[0].parameters);
        free(m.residues);
        assert_int_equal(1, 0);


}

void speedy_gaf(void) {
    clock_t begin = clock();

    //   int GAF_S2(Decimal sig[3], struct Orient ** A, struct Orient ** B, Decimal * S2[], int length,  int mode) {
    int i;
    Decimal S1q = 0, S2q = 0, S3q = 0;
    Decimal *S2[] = {&S1q, &S2q, &S3q};
    Decimal sig[3];
    struct Orient A;
    A.theta = 1.768019;
    A.phi = 3.141593;
    calculate_Y2(&A);

    struct Orient *a[] = {&A, &A, &A};
    Decimal S2sum[3] = {0, 0, 0};
    for (i = 0; i < 10000; i++) {
        sig[0] = (rand() % 100) / 100.;
        sig[1] = (rand() % 100) / 100.;
        sig[2] = (rand() % 100) / 100.;
        GAF_S2(sig, a, a, S2, 3, MODE_REAL);
        S2sum[0] += S1q;
        S2sum[1] += S2q;
        S2sum[2] += S3q;

    }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%lf :: %lf, %lf, %lf\n", time_spent, S2sum[0], S2sum[1], S2sum[2]);


    /* naieve 3.102058 s; */
    /* lookup: 2.141635 s; */

}

int main(void) {
    srand(time(NULL));
    const struct CMUnitTest tests[] = {
            cmocka_unit_test(test_determine_residues),
            cmocka_unit_test(test_temp_tau),
            cmocka_unit_test(test_dwig),
            cmocka_unit_test(test_sphericals),
            cmocka_unit_test(test_rotations),
            cmocka_unit_test(test_gaf),
            cmocka_unit_test(test_statistics),
            cmocka_unit_test(test_crosen_backcalc)
    };
    //speedy_gaf();
    return cmocka_run_group_tests(tests, NULL, NULL);

}
