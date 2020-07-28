#include <stdio.h>
#include "../datatypes.c"
#include "../read_data.c"
#include "../models.c"
//#include "chisq.c"
//#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
//#include "errors.c"
#include <time.h>

long double uniform_rand(void) {
	return ((long double) rand() + 1.) / ((long double) RAND_MAX + 1.);
}


int main(void) {
	srand(time(NULL));
	initialise_dwig();

	char system_file[255] = "system/model.dx";
	int i;
	int err_mod = 0;

	struct Model m;
	int ret = read_system_file(system_file, &m);
	m.error_mode = err_mod;

	if (ret == -1) {
		printf("Error found, crashing peacefully...\n");
		exit(-1);
	}

	
	double sum = 0;
	for (i = 0; i < 10000000; i++) {
		sum += EMF_R1(&(m.residues[40]), &(m.residues[40].relaxation[19]), uniform_rand() * 0.00000001l, uniform_rand(), uniform_rand() * 0.00000000001l, uniform_rand(), MODE_15N);
		sum += EMF_R2(&(m.residues[40]), &(m.residues[40].relaxation[19]), uniform_rand() * 0.00000001l, uniform_rand(), uniform_rand() * 0.00000000001l, uniform_rand(), MODE_15N);
	}
	printf("%f\n", sum);
	
	free_all(&m);
	return 1;
}
