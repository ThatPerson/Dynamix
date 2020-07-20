#include <stdio.h>
#include "datatypes.c"
#include "read_data.c"
#include "models.c"
#include "chisq.c"
#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
#include <time.h>
int main(void) {
	srand(time(NULL));
	struct Model m;
	printf("%d\n", read_system_file("system/model.dx", &m));
	print_system(&m);
	
	int i;
	
	//for (i = 0; i < 1000; i++) {
	for (i = 0; i < m.n_residues; i++) {
		if (m.residues[i].ignore == 1) {
			printf("%d, %f, -1, -1\n", i+1, 1000);
			continue;
		}
		long double opts[2] = {0, 0};
		opts[0] = ((rand() % 100)/100.) * pow(10, -8);
		opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1

		//printf("%Le, %Le\n", opts[0], opts[1]);
		double val = simplex(optimize_chisq, opts, 2, 1.0e-16, 1, &(m.residues[i]), m.model);
		if (val >= 10000) {
			val = -1;
			opts[0] = -1;
			opts[1] = -1;
		}
		printf("%d, %f, %Le, %Le\n", i+1, val, opts[0], opts[1]);
			
	}
	return 1;
}
