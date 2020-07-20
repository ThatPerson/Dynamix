#include <stdio.h>
#include "datatypes.c"
#include "read_data.c"
#include "models.c"
#include "chisq.c"
#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
#include <time.h>
#include <pthread.h>

#define NTHREADS 2


void * run_residue(void *input) {
	int i = ((struct rrargs*)input)->i;
	struct Residue * resid = ((struct rrargs*)input)->resid;
	int model = ((struct rrargs*)input)->model;

	printf("RESIDUE %d\n", i+1);
	printf("Number of relaxations: %f\n", resid->relaxation[1].R);
		
	if (resid->ignore == 1) {
		printf("%d, %f, -1, -1\n", i+1, 1000.);
		return NULL;
	}
	long double opts[2] = {0, 0};
	opts[0] = ((rand() % 100)/100.) * pow(10, -8);
	opts[1] = 0.5 + ((rand() % 100) / 200.); // random number from 0.5 to 1
	//printf("%Le, %Le\n", opts[0], opts[1]);
	double val = simplex(optimize_chisq, opts, 2, 1.0e-16, 1, resid, model);
	if (val >= 10000) {
		val = -1;
		opts[0] = -1;
		opts[1] = -1;
	}
	printf("%d, %f, %Le, %Le\n", i+1, val, opts[0], opts[1]);
}




int main(void) {
	srand(time(NULL));
	struct Model m;
	int ret = read_system_file("system/model.dx", &m);
	printf("%d\n", ret);
	//print_system(&m);
	
	int i;
	
	pthread_t threads[NTHREADS];
	int thread_args[NTHREADS];
	
	int rc;
	
	int current_residue = 0;
	int n_spawns = 0;
	/* if we have 56 residues and 4 threads then we need 
	 * 56 / 4 spawn events (= 14). Add 1 in case (eg for 57).
	 * Then loop over threads, increment current_residue and assign pointers.
	 */
	n_spawns = m.n_residues / NTHREADS;
	n_spawns ++;
	printf("%d spawns\n", n_spawns);

	struct rrargs RRA[NTHREADS]; // = (struct rrargs *)malloc(sizeof(struct rrargs));

	/* spawn the threads */
	int l = 0;
	for (l = 0; l < n_spawns; l++) {
		for (i=0; i<NTHREADS; ++i) {
			RRA[i].i = current_residue + i;
			if (current_residue + i >= m.n_residues)
				continue;
			RRA[i].resid = &(m.residues[current_residue + i]);
			RRA[i].model = m.model;
			thread_args[i] = i;
			printf("spawning thread %d (residue %d)\n", i, current_residue + i);
			rc = pthread_create(&threads[i], NULL, run_residue, (void *) &RRA[i]);
		}
		current_residue += NTHREADS;

	
		for (i=0; i<NTHREADS; ++i) {
			rc = pthread_join(threads[i], NULL);
		}
	}
	//for (i = 0; i < 1000; i++) {
	
	/*for (i = 0; i < m.n_residues; i++) {
		if (m.residues[i].ignore == 1) {
			printf("%d, %f, -1, -1\n", i+1, 1000.);
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
			
	}*/
	free_all(&m);
	return 1;
}
