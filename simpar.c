#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "physics.h"
#include "init_program.h"
#include "debug.h"
#include <time.h>

#include <mpi.h>


/* Fazer um handler para tratar dos argumentos de entrada e passar para o init_particles*/

int main(int argc, char *argv[])
{ 
	int myrank, numberOfProcess;
	MPI_Comm comm;
	MPI_Status status;

	grid_t grid;
	particle_t *par;
	parameters params; 

	MPI_Init( &argc, &argv );
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcess);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// init particles
	params = handler_input(argc, argv);
	
	par = CreateParticleArray(params.n_part);

	if(myrank == 0) {
		init_particles(params, par);
		// Init Total Grid
		grid = initTotalGrid(grid, params.ncside);
		findGridDivision(&params, numberOfProcess);
	}
	else {
		// Init Partial Grid
		grid = initTotalGrid(numberOfProcess, myrank, grid, params.ncside);
	}

	return 0;
}

