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
	int rank, numberOfProcess;
	MPI_Comm comm;
	MPI_Status status;

	grid_t grid;
	grid_tt **gridSendReceive;
	particle_t *par;

	long k;

	MPI_Init( &argc, &argv );
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcess);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// Criação de estrutura particula para enviar em MPI
	const int nitemsPart = 10;
	MPI_Aint displacementsPart[10] = {	offsetof(particle_t, active), 
									offsetof(particle_t, m), 
									offsetof(particle_t, positionX), 
									offsetof(particle_t, positionY),
									offsetof(particle_t, vx), 
									offsetof(particle_t, vy), 
									offsetof(particle_t, gridCoordinateX),
									offsetof(particle_t, gridCoordinateY),
									offsetof(particle_t, appliedForceX), 
									offsetof(particle_t, appliedForceY)};
	int block_lengthsPart[10]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype typesPart[10] = {MPI_INT, MPI_FLOAT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
								MPI_DOUBLE, MPI_LONG, MPI_LONG, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype mpi_particle_t;
	MPI_Type_create_struct(nitemsPart, block_lengthsPart, displacementsPart, typesPart, &mpi_particle_t);
	MPI_Type_commit(&mpi_particle_t);

	// Criação de estrutura grid para enviar em MPI
	const int nitemsGrid = 3;
	MPI_Aint displacementsGrid[3] = {	offsetof(grid_tt, m), 
									offsetof(grid_tt, centerOfMassX), 
									offsetof(grid_tt, centerOfMassY), 
									};
	int block_lengthsGrid[3]  = {1, 1, 1};
	MPI_Datatype typesGrid[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype mpi_grid_t;
	MPI_Type_create_struct(nitemsGrid, block_lengthsGrid, displacementsGrid, typesGrid, &mpi_grid_t);
	MPI_Type_commit(&mpi_grid_t);

	// init particles
	handler_input(argc, argv);
	
	par = CreateParticleArray(params.n_part);

	findGridDivision(numberOfProcess, rank);
	
	grid = initTotalGrid(grid, params.ncside);

	gridSendReceive = initGridSendReceive(rank);
	
	if(rank == 0) {
		// Inicia as particulas e descobre a sua posicao na grelha (pode ser paralelizado a descoberta )
		init_particles(par);
	}
	
	MPI_Bcast(par, params.n_part, mpi_particle_t, 0, MPI_COMM_WORLD);

	for(int i = 0; i < params.n_part; i++) {
		if(par[i].gridCoordinateX >= params.xLowerBound && par[i].gridCoordinateX <= params.xUpperBound && 
		   par[i].gridCoordinateY >= params.yLowerBound && par[i].gridCoordinateY <= params.yUpperBound) {
        	par[i].active = 1;
		}
			
	}


	// Time Step simulation
	for (k = 0; k < params.timeStep; k++) {
		// Run throw all the cells and resets all the center of mass
		memset(grid.m, 0, params.gridSize*sizeof(double));
		memset(grid.centerOfMassX, 0, params.gridSize*sizeof(double));
		memset(grid.centerOfMassY, 0, params.gridSize*sizeof(double));

		for(int i = params.n_part - 1; i >= 0; i--) {
			if(par[i].active == 1) {
				CENTEROFMASSX(par[i].gridCoordinateX, par[i].gridCoordinateY) = CENTEROFMASSX(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m * par[i].positionX;
				CENTEROFMASSY(par[i].gridCoordinateX, par[i].gridCoordinateY) = CENTEROFMASSY(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m * par[i].positionY;
				MASS(par[i].gridCoordinateX, par[i].gridCoordinateY) = MASS(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m;
			}
		}

		int pos[8] = {0};
		for (int i = params.yUpperBound; i >= params.yLowerBound; i--) {
			for (int j = params.xUpperBound; j >= params.xLowerBound; j--) {
				CENTEROFMASSX(i, j) = CENTEROFMASSX(i, j)/MASS(i, j);
				CENTEROFMASSY(i, j) = CENTEROFMASSY(i, j)/MASS(i, j);

				if(j == params.xUpperBound) {
					gridSendReceive[RIGHTPROCESS][pos[RIGHTPROCESS]].centerOfMassX = CENTEROFMASSX(i, j);
					gridSendReceive[RIGHTPROCESS][pos[RIGHTPROCESS]].centerOfMassY = CENTEROFMASSY(i, j);
					gridSendReceive[RIGHTPROCESS][pos[RIGHTPROCESS]++].m = MASS(i, j);
					if(i == params.yUpperBound) {
						gridSendReceive[UPRIGHTPROCESS][pos[UPRIGHTPROCESS]].centerOfMassX = CENTEROFMASSX(i, j);
						gridSendReceive[UPRIGHTPROCESS][pos[UPRIGHTPROCESS]].centerOfMassY = CENTEROFMASSY(i, j);
						gridSendReceive[UPRIGHTPROCESS][pos[UPRIGHTPROCESS]++].m = MASS(i, j);
					}
					if(i == params.yLowerBound)	{
						gridSendReceive[DOWNRIGHTPROCESS][pos[DOWNRIGHTPROCESS]].centerOfMassX = CENTEROFMASSX(i, j);
						gridSendReceive[DOWNRIGHTPROCESS][pos[DOWNRIGHTPROCESS]].centerOfMassY = CENTEROFMASSY(i, j);
						gridSendReceive[DOWNRIGHTPROCESS][pos[DOWNRIGHTPROCESS]++].m = MASS(i, j);
					}
				}

				if(j == params.xLowerBound)	{
					gridSendReceive[LEFTPROCESS][pos[LEFTPROCESS]].centerOfMassX = CENTEROFMASSX(i, j);
					gridSendReceive[LEFTPROCESS][pos[LEFTPROCESS]].centerOfMassY = CENTEROFMASSY(i, j);
					gridSendReceive[LEFTPROCESS][pos[LEFTPROCESS]++].m = MASS(i, j);
					if(i == params.yUpperBound) {
						gridSendReceive[UPLEFTPROCESS][pos[UPLEFTPROCESS]].centerOfMassX = CENTEROFMASSX(i, j);
						gridSendReceive[UPLEFTPROCESS][pos[UPLEFTPROCESS]].centerOfMassY = CENTEROFMASSY(i, j);
						gridSendReceive[UPLEFTPROCESS][pos[UPLEFTPROCESS]++].m = MASS(i, j);
					}
					if(i == params.yLowerBound)	{
						gridSendReceive[DOWNLEFTPROCESS][pos[DOWNLEFTPROCESS]].centerOfMassX = CENTEROFMASSX(i, j);
						gridSendReceive[DOWNLEFTPROCESS][pos[DOWNLEFTPROCESS]].centerOfMassY = CENTEROFMASSY(i, j);
						gridSendReceive[DOWNLEFTPROCESS][pos[DOWNLEFTPROCESS]++].m = MASS(i, j);
					}
				}

				if(i == params.yUpperBound) {
					gridSendReceive[UPPROCESS][pos[UPPROCESS]].centerOfMassX = CENTEROFMASSX(i, j);
					gridSendReceive[UPPROCESS][pos[UPPROCESS]].centerOfMassY = CENTEROFMASSY(i, j);
					gridSendReceive[UPPROCESS][pos[UPPROCESS]++].m = MASS(i, j);
				}

				if(i == params.yLowerBound)	{
					gridSendReceive[RIGHTPROCESS][pos[RIGHTPROCESS]].centerOfMassX = CENTEROFMASSX(i, j);
					gridSendReceive[RIGHTPROCESS][pos[RIGHTPROCESS]].centerOfMassY = CENTEROFMASSY(i, j);
					gridSendReceive[RIGHTPROCESS][pos[RIGHTPROCESS]++].m = MASS(i, j);
				}
			}
		}

		// Organizado em grelha
		if(params.xSize != 1) {
			// Enviar à esquerda
			if(params.xLowerBound == 0)
				MPI_Irecv(gridSendReceive[8], params.sizeSmallD, mpi_grid_t, rank + params.xSize - 1, MYTAG, WORLD, &request);
			else
				MPI_Irecv(gridSendReceive[8], params.sizeSmallD, mpi_grid_t, rank - 1, MYTAG, WORLD, &request);

				// Enviar à esquerda cima
				if(params.xLowerBound == 0)
					MPI_Irecv(gridSendReceive[9], params.sizeSmallD, mpi_grid_t, rank + 2*params.xSize - 1, MYTAG, WORLD, &request);
				else
					MPI_Irecv(gridSendReceive[9], params.sizeSmallD, mpi_grid_t, rank + params.xSize, MYTAG, WORLD, &request);

			// Enviar à cima
			if(params.yUpperBound == params.ncside - 1)
				MPI_Irecv(gridSendReceive[10], params.sizeSmallD, mpi_grid_t, rank - numberOfProcess + params.xSize, MYTAG, WORLD, &request);
			else
				MPI_Irecv(gridSendReceive[10], params.sizeSmallD, mpi_grid_t, rank + params.xSize, MYTAG, WORLD, &request);

				// Enviar à direita cima
				if(params.xLowerBound == 0)
					MPI_Irecv(gridSendReceive[11], params.sizeSmallD, mpi_grid_t, rank + params.xSize - 1, MYTAG, WORLD, &request);
				else
					MPI_Irecv(gridSendReceive[11], params.sizeSmallD, mpi_grid_t, rank - 1, MYTAG, WORLD, &request);

			// Enviar à direita
			if(params.xUpperBound == params.ncside - 1)
				MPI_Irecv(gridSendReceive[12], params.sizeSmallD, mpi_grid_t, rank - params.xSize + 1, MYTAG, WORLD, &request);
			else
				MPI_Irecv(gridSendReceive[12], params.sizeSmallD, mpi_grid_t, rank + 1, MYTAG, WORLD, &request);

				// Enviar à direita baixo
				if(params.xLowerBound == 0)
					MPI_Irecv(gridSendReceive[13], params.sizeSmallD, mpi_grid_t, rank + 2*params.xSize - 1, MYTAG, WORLD, &request);
				else
					MPI_Irecv(gridSendReceive[13], params.sizeSmallD, mpi_grid_t, rank + params.xSize, MYTAG, WORLD, &request);

			// Enviar à baixo
			if(params.yLowerBound == 0)
				MPI_Irecv(gridSendReceive[14], params.sizeSmallD, mpi_grid_t, rank + numberOfProcess - params.xSize, MYTAG, WORLD, &request);
			else
				MPI_Irecv(gridSendReceive[14], params.sizeSmallD, mpi_grid_t, rank - params.xSize, MYTAG, WORLD, &request);

				// Enviar à esquerda baixo
				if(params.xLowerBound == 0)
					MPI_Irecv(gridSendReceive[15], params.sizeSmallD, mpi_grid_t, rank + params.xSize - 1, MYTAG, WORLD, &request);
				else
					MPI_Irecv(gridSendReceive[15], params.sizeSmallD, mpi_grid_t, rank - 1, MYTAG, WORLD, &request);


			MPI_Send(A, 100, MPI_DOUBLE, 1, MYTAG, WORLD);
			MPI_Wait(&request, &status);
		}
		// Organizada em linhas
		else {

		}
/*
		// Compute interactions
		// Run all particles
		for(int i = params.n_part - 1; i >= 0; i = i - 1){
			if(par[i].active == 1) {
				par[i].appliedForceX = 0;
				par[i].appliedForceY = 0;

				// Run the adjacent grids
				int j;
				int m;
				for (int n = 0; n < 9; n = n + 1) {
					j = n / 3 - 1;
					m = n % 3 - 1;
					aux1 = 0, aux2 = 0;

					if(par[i].gridCoordinateX+j == -1) {
						sideLEFTRIGHT = LEFT;
						aux1 = par[i].gridCoordinateX+j + params.ncside;
					} else if(par[i].gridCoordinateX+j == params.ncside) {
						sideLEFTRIGHT = RIGHT;
						aux1 = par[i].gridCoordinateX+j - params.ncside;
					} else {
						sideLEFTRIGHT = MIDDLE;
						aux1 = par[i].gridCoordinateX+j;
					}

					if(par[i].gridCoordinateY+m == -1){
						sideUPDOWN = DOWN;
						aux2 = par[i].gridCoordinateY+m + params.ncside;
					} else if(par[i].gridCoordinateY+m == params.ncside){
						sideUPDOWN = UP;
						aux2 = par[i].gridCoordinateY+m - params.ncside;
					} else {
						sideUPDOWN = MIDDLE;
						aux2 = par[i].gridCoordinateY+m;
					}
					
					
					if(grid.m[ MATRIX(aux1, aux2, params.ncside) ] != 0)
						calculateGravForce(&(par[i]), CENTEROFMASSX(aux1, aux2), CENTEROFMASSY(aux1, aux2), 
												MASS(aux1, aux2), sideUPDOWN, sideLEFTRIGHT); //for each adjacent cell.---- 
				}


				invM = 1.0/par[i].m;
				// Updates particles position and velocity and position on the grid
				par[i].vx += par[i].appliedForceX * invM; //a = F/m 
				par[i].vy += par[i].appliedForceY * invM;

				par[i].positionX += par[i].vx + 0.5 * par[i].appliedForceX * invM;//x = x0 + v0t + 0.5 a t^2 (t = 1)
				par[i].positionY += par[i].vy + 0.5 * par[i].appliedForceY * invM;

				//See if its out of bounds
				if(par[i].positionX >= 1) par[i].positionX = par[i].positionX - (int)(par[i].positionX);
				else if(par[i].positionX < 0) par[i].positionX = 1 + (par[i].positionX - ceil(par[i].positionX)); 

				if(par[i].positionY >= 1) par[i].positionY = par[i].positionY - (int)(par[i].positionY);
				else if(par[i].positionY < 0) par[i].positionY = 1 + (par[i].positionY - ceil(par[i].positionY));

				// Updates the position of the particle on the grid of cells
				par[i].gridCoordinateX = par[i].positionX * params.ncside;
				par[i].gridCoordinateY = par[i].positionY * params.ncside;
		}
			}*/
	}

	MPI_Type_free(&mpi_particle_t);
	MPI_Finalize();

	return 0;
}

