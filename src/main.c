/**
 * @file main.c
 */

#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "datatypes.h"
#include "read_data.h"

#include "models/smf.h"
#include "models/emf.h"
#include "models/gaf.h"
#include "models/egaf.h"
#include "models/aimf.h"
#include "chisq.h"
#include "crosen.h" // implementation of Nelder-Mead simplex algorithm
#include "anneal.h" // simulated annealing algorithm
#include "errors.h"
#include "global_gaf.h"
#include "runners.h"





#include "verification.h"

#include "runners.h"



/**
 * Start function. Initialises parameters, loads files, spawns threads and outputs data.
 * @opts filename
 *  Takes input as path to file containing System definition (*dx file)
 * @opts -e
 *  Enables error mode
 */
int main(int argc, char * argv[]) {
	MPI_Init( &argc, &argv );
	int numprocs, myid;
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs ); 
	MPI_Comm_rank( MPI_COMM_WORLD, &myid );
	
	//printf("I am %d/%d.\n", myid, numprocs);

	/* Initialisation */
	start_time = time(0);
	srand((unsigned int)time(NULL));
	initialise_dwig(HALF_PI, Dwig);


	char system_file[255] = "";
	unsigned int i;
	int err_mod = 0;
	for (i = 1; i < (unsigned int) argc; i++) {
		//printf("%s\n", argv[i]);
		if (strcmp(argv[i], "-e") == 0)
			err_mod = 1;
		else if (strcmp(argv[i], "verify") == 0)
			verify_all();
		else
			strcpy(system_file, argv[i]);
	}

	if (strcmp(system_file, "") == 0){
		printf("Please provide system file.\n");
		exit(-1);
	}



	struct Model m;
	int ret = read_system_file(system_file, &m);
    if (ret == -1) {
        printf("Error found, crashing peacefully...\n");
        free_all(&m);
        MPI_Finalize();
        exit(-1);
    }
	m.myid = myid;
	m.numprocs = numprocs;
	//printf("%d params.\n", m.params);
	//exit(-1);
	determine_residues(&m, myid, numprocs, &(m.proc_start), &(m.proc_end));
	printf("%d/%d: running %d - %d\n", myid+1, numprocs, m.proc_start, m.proc_end);
	omp_set_num_threads(m.nthreads);
	
	m.error_mode = err_mod;
	if (m.error_mode == 1 && m.n_error_iter == 0) {
		printf("Please provide number of error iterations\n");
		ret = -1;
	}

	if (m.model == MOD_UNDEFINED) {
		printf("Model undefined\n");
		ret = -1;
	}

    if (ret == -1) {
        free_all(&m);
        MPI_Finalize();
    }


	//printf("%d\n", ret);
	if (m.myid == 0) {
		char filename[300];
		sprintf(filename, "%s/model.txt", m.outputdir);
		print_system(&m, filename);
	}
	if (m.global == LOCAL)
		run_fitting(&m);
	else if (m.global == GLOBAL)
		run_global(&m);
	
	
	
	printf("Worker %d: Printing residues...\n", m.myid+1);
	print_residues(&m);
	
	if (m.error_mode == 1) {
		printf("Worker %d: Running errors...\n", m.myid+1);
		if (m.global == LOCAL)
			run_errors(&m);
		else
			calc_global_errors(&m);
		printf("Worker %d: Printing errors...\n", m.myid+1);
		print_errors(&m);
	}
	
	printf("Worker %d: Printing gaf...\n", m.myid+1);
	print_gaf(&m);
	print_backcalcs(&m);
	


	free_all(&m);
	MPI_Finalize();
	return 0;
}
