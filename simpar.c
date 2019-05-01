#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "physics.h"
#include "init_program.h"
#include "debug.h"
#include <time.h>
#include <stddef.h>

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

	long k;

	MPI_Init( &argc, &argv );
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcess);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Criação de estrutura particula para enviar em MPI
	MPI_Datatype particle_dt;
	MPI_Aint displacements[10] = {	offsetof(particle_t, active), 
									offsetof(particle_t, m), 
									offsetof(particle_t, positionX), 
									offsetof(particle_t, positionY),
									offsetof(particle_t, vx), 
									offsetof(particle_t, vy), 
									offsetof(particle_t, gridCoordinateX),
									offsetof(particle_t, gridCoordinateY),
									offsetof(particle_t, appliedForceX), 
									offsetof(particle_t, appliedForceY)};
	int block_lengths[10]  = {1};
	MPI_Datatype types[10] = {MPI_INT, MPI_FLOAT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
								MPI_DOUBLE, MPI_LONG, MPI_LONG, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Type_create_struct(10, block_lengths, displacements, types, &particle_dt);
	MPI_Type_commit(&particle_dt);

	// init particles
	params = handler_input(argc, argv);
	
	par = CreateParticleArray(params.n_part);

	findGridDivision(&params, numberOfProcess);
	
	grid = initTotalGrid(grid, params.ncside);
	
	if(myrank == 0) {
		// Inicia as particulas e descobre a sua posicao na grelha (pode ser paralelizado a descoberta )
		init_particles(params, par);
		
		// Divisão e envio das particulas pelos varios processos
	}
	//else {
		// Init Partial Grid
		//grid = initPartialGrid(numberOfProcess, myrank, grid, &params);

		// Recessão das particulas pelos varios processos
			// Contando o numero de particulas recebidas para atualziar params.partialNrPart
	//}

	MPI_Bcast(par, params.n_part*sizeof(particle_t), particle_dt, 0, MPI_COMM_WORLD);


	// Time Step simulation
	for (k = 0; k < params.timeStep; k++) {
		// Run throw all the cells and resets all the center of mass
		memset(grid.m, 0, params.gridSize*sizeof(double));
		memset(grid.centerOfMassX, 0, params.gridSize*sizeof(double));
		memset(grid.centerOfMassY, 0, params.gridSize*sizeof(double));

		//long x, y;

		for(int i = 0; i < params.partialNrPart; i++) {
			//x = par[i].gridCoordinateX;
			//y = par[i].gridCoordinateY;
		}
	}

	return 0;
}

