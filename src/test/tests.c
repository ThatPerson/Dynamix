//
// Created by ben on 08/03/2021.
//

#include "tests.h"
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

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
    assert_int_equal(start, 42);
    assert_int_equal(end,56 );
    k = determine_residues(&m, 2, 3, &start, &end);
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

int main(void) {
    const struct CMUnitTest tests[] = {
            cmocka_unit_test(test_determine_residues),
            cmocka_unit_test(test_temp_tau),
            cmocka_unit_test(test_dwig),
            cmocka_unit_test(test_sphericals),
            cmocka_unit_test(test_rotations),
            cmocka_unit_test(test_gaf),
            cmocka_unit_test(test_statistics)
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
