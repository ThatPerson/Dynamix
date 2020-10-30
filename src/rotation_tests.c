#include <stdio.h>
#include "datatypes.h"
#include <math.h>

#include <time.h>



int main(int argc, char * argv[]) {

	/* Initialisation */
	srand((unsigned int)time(NULL));
	initialise_dwig(HALF_PI, Dwig);

	int ml, mi;
	/*for (ml = 2; ml >= -2; ml--) {
		for (mi = 2; mi >= -2; mi--) {
			printf("%Lf\t\t", Dwig[mi+2][ml+2], mi, ml);
			// m' is mi+2
			// m  is ml+2

		}
		printf("\n");
	}*/

	int m, mp;
	printf("Rows are m', columns are m\n");
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


	struct Orient test;
	test.theta = 0.1;
	test.phi = 0.2;
	calculate_Y2(&test);
	for (mi = -2; mi <= 2; mi++) {
		printf("%d: %lf + %lf I\n", mi, creal(test.Y2[mi+2]), cimag(test.Y2[mi+2]));
	}

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
	int i;
	printf("CACA: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(CACA.Y2[i]), cimag(CACA.Y2[i]));
	}
	printf("pp : \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(pp.Y2[i]), cimag(pp.Y2[i]));
	}
	printf("NHv: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(NHv.Y2[i]), cimag(NHv.Y2[i]));
	}
	//void rotate_Y2(struct Orient * or, double alpha, double beta, double gamma) {
	printf("Rotating...\n");
	rotate_Y2(&CACA, HALF_PI, M_PI, HALF_PI);
	rotate_Y2(&pp, HALF_PI, M_PI, HALF_PI);
	rotate_Y2(&NHv, HALF_PI, M_PI, HALF_PI);
	printf("CACA: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(CACA.Y2[i]), cimag(CACA.Y2[i]));
	}
	printf("pp : \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(pp.Y2[i]), cimag(pp.Y2[i]));
	}
	printf("NHv: \n");
	for (i = 0; i < 5; i++) {
		printf("\t%d:\t%lf + %lf I\n", i-2, creal(NHv.Y2[i]), cimag(NHv.Y2[i]));
	}
	return 1;
}
