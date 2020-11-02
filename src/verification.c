/**
 * @file verification.c
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "datatypes.h"

//__func__,__FILE__, __LINE__

#define LOGGING

#ifdef LOGGING
#define error(a, args...) printf("ERROR %s(%s:%d %d): " a "\n",  __func__,__FILE__, __LINE__, (int) (time(0) - start_time), ##args)
#define warn(a, args...) printf("WARN  %s(%s:%d %d): " a "\n",  __func__,__FILE__, __LINE__, (int) (time(0) - start_time), ##args)
#define log(a, args...) printf("LOG   %s(%s:%d %d): " a "\n",  __func__,__FILE__, __LINE__, (int) (time(0) - start_time), ##args)
#else
#define error(a, args...) 
#define warn(a, args...)
#define log(a, args...)
#endif

#define N_BINS 10

void verify_rotations(void) {
	printf("\n\n\n===== VERIFYING ROTATION =====\n");
	printf(">>> \t*Verifying Wigner D Matrix*\n");
	printf(">>> \tI will output the pi/2 reduced Wigner D matrix.\n");

	int m, mp, i;
	printf("Rows are m', columns are m\n\n");
	printf("\t\t");
	for (m = 2; m >= -2; m--) {
		printf("%d\t", m);
	}
	printf("\n");
	for (mp = 2; mp >= -2; mp--) {
		printf("%d\t", mp);
		for (m = 2; m >= -2; m--) {
			printf("%Lf\t", Dwig[mp+2][m+2]);
		}
		printf("\n");
	}

	printf("\n\n");

	printf(">>> \tThe following Mathematica code should output the same.\n\n");
	printf("d[j_, mb_, m_, beta_] := \nSqrt[Factorial [j + mb] Factorial[j - mb] Factorial[j + m] Factorial[\n     j - m]] Sum[((-1)^(mb - m + s) (Cos[beta/2])^(\n    2 j + m - mb - 2 s) (Sin[beta/2])^(mb - m + 2 s))/(\n   Factorial[j + m - s] Factorial[s] Factorial[mb - m + s] Factorial[\n     j - mb - s]), {s, 0, j - mb}]\nTable[{mp, m, N[ d[2, mp, m, Pi/2]]}, {mp, 2, -2, -1}, {m, \n   2, -2, -1}] // MatrixForm\n\n");
	printf(">>> \tThis reduced Wigner function comes from:\n\t\thttps://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_(small)_d-matrix\n");

	float beta = (rand() % 100) / 100.;

	printf(">>> \tI will now output the reduced Wigner D matrix for beta = %f. \n", beta);

	long double test_dwig[5][5];
	initialise_dwig(beta, test_dwig);

	printf(">>> Rows are m', columns are m\n\n");
	printf("\t\t");
	for (m = 2; m >= -2; m--) {
		printf("%d\t", m);
	}
	printf("\n");
	for (mp = 2; mp >= -2; mp--) {
		printf("%d\t", mp);
		for (m = 2; m >= -2; m--) {
			printf("%Lf\t", test_dwig[mp+2][m+2]);
		}
		printf("\n");
	}

	printf("\n\n");

	printf(">>> \tThe following Mathematica code should output the same.\n\n");
	printf("d[j_, mb_, m_, beta_] := \nSqrt[Factorial [j + mb] Factorial[j - mb] Factorial[j + m] Factorial[\n     j - m]] Sum[((-1)^(mb - m + s) (Cos[beta/2])^(\n    2 j + m - mb - 2 s) (Sin[beta/2])^(mb - m + 2 s))/(\n   Factorial[j + m - s] Factorial[s] Factorial[mb - m + s] Factorial[\n     j - mb - s]), {s, 0, j - mb}]\nTable[{mp, m, N[ d[2, mp, m, %f]]}, {mp, 2, -2, -1}, {m, \n   2, -2, -1}] // MatrixForm\n\n", beta);
	
	printf(">>> *Verifying Spherical Harmonics.*\n");

	struct Orient test;
	test.theta = (rand() % 628)/100.;
	test.phi = (rand() % 628)/100.;
	printf(">>> \tI have calculated the Y2m set for theta = %f, phi = %f. \n\n", test.theta, test.phi);
	calculate_Y2(&test);
	for (m = -2; m <= 2; m++) {
		printf("Y^2_%d:  %lf + %lf I\n", m, creal(test.Y2[m+2]), cimag(test.Y2[m+2]));
	}
	printf(">>> \tAnd conjugates...\n");
	for (m = -2; m <= 2; m++) {
		printf("Y^c2_%d: %lf + %lf I\n", m, creal(test.Y2c[m+2]), cimag(test.Y2c[m+2]));
	}

	printf("\n>>> \tThis Mathematica code should give the same.\n\n");
	printf("Table[{m, SphericalHarmonicY[2, m, %f, %f], \n   Conjugate[SphericalHarmonicY[2, m, %f, %f]]}, {m, -2, \n   2}] // TableForm\n\n", test.theta, test.phi, test.theta, test.phi);

	printf(">>> \tPlease note that the original MATLAB version had the spherical harmonics the wrong way around (though in practice this made little difference).\n");

	printf(">>> *Verifying Rotations of SH.*\n");
	printf(">>> \tI will create three orientations; \n");

	struct Orient CACA; // ca ca axis
	struct Orient pp; // plane perpendicula
	struct Orient NHv; // roughly parallel to N-H

	CACA.theta = 0;
	CACA.phi = 0;
	pp.theta = HALF_PI;
	pp.phi = HALF_PI;
	NHv.theta = HALF_PI ;
	NHv.phi = 0;

	calculate_Y2(&CACA);
	calculate_Y2(&pp);
	calculate_Y2(&NHv);

	printf(">>> \t\tA, oriented along the Z (gamma) axis (%f, %f)\n", CACA.theta, CACA.phi);
	printf(">>> \t\tB, oriented along the Y (beta) axis (%f, %f)\n", pp.theta, pp.phi);
	printf(">>> \t\tC, oriented along the X (alpha) (%f, %f)\n", NHv.theta, NHv.phi);
	printf(">>> \tThese have the following Spherical Harmonic components;\n\n");

	printf("A: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(CACA.Y2[i]), cimag(CACA.Y2[i]));
	}
	printf("B : \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(pp.Y2[i]), cimag(pp.Y2[i]));
	}
	printf("C: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(NHv.Y2[i]), cimag(NHv.Y2[i]));
	}
	//void rotate_Y2(struct Orient * or, double alpha, double beta, double gamma) {

	printf("\n>>> \tApplying a rotation of pi/2 about the Z axis should affect B & C but not A\n\n");
	rotate_Y2(&CACA, HALF_PI, 0, 0);
	rotate_Y2(&pp, HALF_PI, 0, 0);
	rotate_Y2(&NHv, HALF_PI, 0, 0);
	printf("A: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(CACA.Y2[i]), cimag(CACA.Y2[i]));
	}
	printf("B : \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(pp.Y2[i]), cimag(pp.Y2[i]));
	}
	printf("C: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(NHv.Y2[i]), cimag(NHv.Y2[i]));
	}

	calculate_Y2(&CACA);
	calculate_Y2(&pp);
	calculate_Y2(&NHv);

	printf("\n>>> \tApplying a rotation of pi/2 about the Y axis should affect A & C but not B\n\n");

	rotate_Y2(&CACA, 0, HALF_PI, 0);
	rotate_Y2(&pp, 0, HALF_PI, 0);
	rotate_Y2(&NHv, 0, HALF_PI, 0);
	printf("A: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(CACA.Y2[i]), cimag(CACA.Y2[i]));
	}
	printf("B : \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(pp.Y2[i]), cimag(pp.Y2[i]));
	}
	printf("C: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(NHv.Y2[i]), cimag(NHv.Y2[i]));
	}

	calculate_Y2(&CACA);
	calculate_Y2(&pp);
	calculate_Y2(&NHv);

	printf("\n>>> \tApplying a rotation of pi/2 about Y followed by pi/2 about Z should affect all, but A should be as after pi/2 about Y\n\n");

	rotate_Y2(&CACA, 0, HALF_PI, HALF_PI);
	rotate_Y2(&pp, 0, HALF_PI, HALF_PI);
	rotate_Y2(&NHv, 0, HALF_PI, HALF_PI);
	printf("A: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(CACA.Y2[i]), cimag(CACA.Y2[i]));
	}
	printf("B : \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(pp.Y2[i]), cimag(pp.Y2[i]));
	}
	printf("C: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(NHv.Y2[i]), cimag(NHv.Y2[i]));
	}

	error("Test");

	return;
}

void verify_GAF(void) {
	printf("\n\n\n===== VERIFYING GAF =====\n");
	long double sigs[3];
	int i;
	for (i = 0; i < 3; i++) {
		 sigs[i] = (rand()%100)/100.;
	}
	printf(">>> \tTaking sA = %f, sB = %f, sG = %f\n", (double) sigs[0], (double) sigs[1], (double) sigs[2]);
	struct Orient A, B;
	A.theta = (rand() % 628)/100.;
	A.phi = (rand() % 314)/100.;
	B.theta = (rand() % 628)/100.;
	B.phi = (rand() % 314)/100.;
	calculate_Y2(&A);
	calculate_Y2(&B);
	printf(">>> \tCreating two vectors, A (%f, %f) and B (%f, %f)", A.theta, A.phi, B.theta, B.phi);
	double SAA, SBB, SAB;
	printf(">>> \tI will now calculate real mode order parameters for AxA, BxB, and AxB\n\n");
	struct Orient *set1[] = {&A, &B, &A};
	struct Orient *set2[] = {&A, &B, &B};
	double *S2[] = {&SAA, &SBB, &SAB};

	GAF_S2(sigs, set1, set2, S2, 3, MODE_REAL);

	printf("AxA = %f\nBxB = %f\nAxB = %f\n\n", SAA, SBB, SAB);

	printf(">>> \tI will now output a python file, verify.py, which should give the same values.\n");
	FILE * fp;
	fp = fopen("verify.py", "w");
	fprintf(fp, "import numpy as np\n\ndef sdwig22(beta):\n\tdwig = np.zeros((5, 5))\n\tcb = np.cos(beta)\n\tsb = np.sin(beta)\n\tdwig[0, 0] = np.power(np.cos(beta/2.), 4.)\n\tdwig[1, 0] = (-1/2.) * (1 + cb) * sb\n\tdwig[2, 0] = np.sqrt(3/8.) * np.power(sb, 2.)\n\tdwig[3, 0] = (1/2.) * sb * (cb - 1)\n\tdwig[4, 0] = np.power(np.sin(beta / 2.), 4.)\n\n\tdwig[0, 1] = (1/2.) * (1 + cb) * sb\n\tdwig[1, 1] = np.power(cb, 2.) - (1 - cb)/2.\n\tdwig[2, 1] = -np.sqrt(3/8.) * np.sin(2*beta)\n\tdwig[3, 1] = (1/2.) * (1 + cb) - (cb * cb)\n\tdwig[4, 1] = (1/2.) * (cb - 1) * sb\n\n\tdwig[0, 2] = np.sqrt(3/8.) * (sb * sb)\n\tdwig[1, 2] = np.sqrt(3/2.) * sb * cb\n\tdwig[2, 2] = (1/2. ) * (3 * cb * cb - 1)\n\tdwig[3, 2] = -np.sqrt(3/2.) * sb * cb\n\tdwig[4, 2] = np.sqrt(3/8.) * sb * sb\n\n\tdwig[0, 3] = (-1/2.) * (cb - 1) * sb\n\tdwig[1, 3] = (1/2.) * (1 + cb) - (cb * cb)\n\tdwig[2, 3] = np.sqrt(3/8.) * np.sin(2 * beta)\n\tdwig[3, 3] = (cb * cb) - (1 - cb)/2.\n\tdwig[4, 3] = (-1/2.) * (1 + cb) * sb\n\n\tdwig[0, 4] = np.power(np.sin(beta/2.), 4.)\n\tdwig[1, 4] = (-1/2.) * (cb - 1) * sb\n\tdwig[2, 4] = np.sqrt(3/8.) * sb * sb\n\tdwig[3, 4] = (1/2.) * sb * (cb + 1)\n\tdwig[4, 4] = np.power(np.cos(beta/2.), 4)\n\n\treturn dwig\n\ndwig = sdwig22(np.pi/2.)\ncount = 0\n\ndef Y2mm(theta, phi):\n\tY = np.zeros(5, dtype=complex)\n\tst = np.sin(theta)\n\tct = np.cos(theta)\n\n\tY[0] = (1/4.) * np.sqrt(15 / (2 * np.pi)) * (st * st) * np.exp(-2j * phi)\n\tY[1] = (1/2.) * np.sqrt(15 / (2 * np.pi)) * st * ct * np.exp(-1j * phi)\n\tY[2] = (1/4.) * np.sqrt(5./np.pi) * (3 * ct * ct - 1)\n\tY[3] = (-1/2.) * np.sqrt(15 / (2 * np.pi)) * st * ct * np.exp(1j * phi)\n\tY[4] = (1/4.) * np.sqrt(15 / (2 * np.pi)) * st * st * np.exp(2j * phi)\n\n\treturn Y\n\n\ndef GAF_ord_paramTT(sigA, sigB, sigG, cphi1, ctheta1, cphi2, ctheta2):\n\tglobal dwig\n\tglobal count\n\tcount = count + 1\n\tYY = Y2mm(ctheta1, cphi1)\n\tYYc = np.conj(Y2mm(ctheta2, cphi2))\n\tA = 0\n\tfor l in range(-2, 3):\n\t\tfor m in range(-2, 3):\n\t\t\tfor mp in range(-2, 3): \n\t\t\t\tfor k in range(-2, 3):\n\t\t\t\t\tfor kp in range(-2, 3):\n\t\t\t\t\t\ttemp = np.power(-1j, k - kp)\n\t\t\t\t\t\texps = sigA * sigA * (k**2 + kp**2) / 2.\n\t\t\t\t\t\texps += sigB * sigB * l * l\n\t\t\t\t\t\texps += sigG * sigG * (m**2 + mp**2) / 2.\n\t\t\t\t\t\ttemp *= np.exp(-exps)\n\t\t\t\t\t\ttemp *= dwig[k + 2, l + 2]\n\t\t\t\t\t\ttemp *= dwig[kp + 2, l + 2]\n\t\t\t\t\t\ttemp *= dwig[m + 2, k + 2]\n\t\t\t\t\t\ttemp *= dwig[mp + 2, kp + 2]\n\t\t\t\t\t\ttemp *= YY[m+2]\n\t\t\t\t\t\ttemp *= YYc[mp + 2]\n\t\t\t\t\t\tA = A + temp\n\treturn np.real(A*(4 * np.pi)/5.)\n\nAtheta = %f\nAphi = %f\nBtheta = %f\nBphi = %f\nsA = %f\nsB = %f\nsG = %f\n\nprint('AxA = %%f' %% (GAF_ord_paramTT(sA, sB, sG, Aphi, Atheta, Aphi, Atheta)))\nprint('BxB = %%f' %% (GAF_ord_paramTT(sA, sB, sG, Bphi, Btheta, Bphi, Btheta)))\nprint('AxB = %%f' %% (GAF_ord_paramTT(sA, sB, sG, Aphi, Atheta, Bphi, Btheta)))\n\n\n\n\n\n", A.theta, A.phi, B.theta, B.phi, (double) sigs[0], (double) sigs[1], (double) sigs[2]);

	fclose(fp);

	return;
}

void verify_stats(void) {
	printf("\n\n\n===== VERIFYING STATS =====\n");
	int i;
	int bins[N_BINS];
	for (i = 0; i < N_BINS; i++) {
		bins[i] = 0;
	}
	int val;
	printf(">>> \tI will generate a sample of 40000 uniform random numbers between 0 and 1. I will place these in %d bins.\n\n", N_BINS);
	for (i = 0; i < 40000; i++) {
		val = (int) (uniform_rand() * N_BINS);
		bins[val]++;
	}
	for (i = 0; i < N_BINS; i++) {
		printf("Bin %d: %d\n", i, bins[i]);
	}
	printf(">>> \tPlease verify these appear uniform.\n");
	double mean = (rand()%1000)/1000.;
	double std = (rand()%1000)/1000.;
	long double calc_mean, calc_std;
	printf(">>> \tI will now generate a list of 20000 normally distributed values with mean %f and std %f using norm_rand()\n", mean, std);
	printf(">>> \tI will then calculate statistics using calc_statistics.\n");

	long double rnds[20000];
	for (i = 0; i < 20000; i++) {
		rnds[i] = (long double) norm_rand(mean, std);
	}
	calc_statistics(rnds, 20000, &calc_mean, &calc_std);
	printf(">>> \tCalculated mean: %Lf, std: %Lf\n", calc_mean, calc_std);

	//long double uniform_rand(void);
	//double norm_rand(double mean, double std);
	//void calc_statistics(long double * vals, int length, long double * mean, long double * std);
}

void verify_all(void) {
	verify_rotations();
	verify_GAF();
	verify_stats();
	return;
}