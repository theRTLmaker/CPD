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

#define SENDCENTER 0

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
	int idToSend[8] = {0};
	int flag;

	MPI_Request request;

	MPI_Init( &argc, &argv );
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcess);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// Criação de estrutura particula completa para enviar em MPI
	const int nitemsPart = 10;
	MPI_Aint displacementsPart[10] = {	offsetof(particle_t, number), 
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
	MPI_Datatype typesPart[10] = {MPI_LONG_LONG_INT, MPI_FLOAT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
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
        	par[i].number = i;
		}
	}		

/*  						topo
	lateral esquerdo	-------------  lateral direito
				    	| 1 | 2 | 3 |
					    -------------
						| 0 |   | 4 |
						-------------
						| 7 | 6 | 5 |
						-------------
							baixo

*/
	// Organizado em grelha
	if(params.xSize != 1) {

		// Enviar à esquerda
		if(params.xLowerBound == 0)
			idToSend[0] = rank + params.xSize - 1;
		else
			idToSend[0] = rank - 1;

		// Enviar à esquerda cima
		// Canto superior esquerdo
		if(params.xLowerBound == 0 && params.yUpperBound == params.ncside - 1)
			idToSend[1] = params.xSize - 1;
		// Lateral esquerdo
		else if(params.xLowerBound == 0 && params.yUpperBound != params.ncside - 1)
			idToSend[1] = rank + 2*params.xSize - 1;
		// Topo
		else if(params.xLowerBound != 0 && params.yUpperBound == params.ncside - 1)
			idToSend[1] = rank - numberOfProcess + params.xSize - 1;
		else
			idToSend[1] = rank + params.xSize - 1;

		// Enviar à cima
		if(params.yUpperBound == params.ncside - 1)
			idToSend[2] = rank - numberOfProcess + params.xSize;
		else
			idToSend[2] =  rank + params.xSize;

		// Enviar à direita cima
		// Canto superior direito
		if (params.xUpperBound == params.ncside - 1 && params.yUpperBound == params.ncside - 1)
			idToSend[3] =  0;
		// Lateral direito
		else if(params.xUpperBound == params.ncside - 1 && params.yUpperBound != params.ncside - 1)
			idToSend[3] =  rank + 1;
		// Topo
		else if(params.xUpperBound != params.ncside - 1 && params.yUpperBound == params.ncside - 1)
			idToSend[3] =  rank - numberOfProcess + params.xSize + 1;
		else
			idToSend[3] =  rank + params.xSize + 1;

		// Enviar à direita
		if(params.xUpperBound == params.ncside - 1)
			idToSend[4] = rank - params.xSize + 1;
		else
			idToSend[4] = rank + 1;

		// Enviar à direita baixo
		// Canto inferior direito
		if(params.xUpperBound == params.ncside - 1 && params.yLowerBound == 0)
			idToSend[5] = numberOfProcess - params.xSize;
		// Baixo
		else if(params.xUpperBound != params.ncside - 1 && params.yLowerBound == 0)
			idToSend[5] =  rank + numberOfProcess - params.xSize + 1;
		// Lateral direito
		else if(params.xUpperBound == params.ncside - 1 && params.yLowerBound != 0)
			idToSend[5] =  rank - 2*params.xSize + 1;
		else
			idToSend[5] = rank - params.xSize + 1;

		// Enviar a baixo
		if(params.yLowerBound == 0)
			idToSend[6] = rank + numberOfProcess - params.xSize;
		else
			idToSend[6] = rank - params.xSize;

		// Enviar à esquerda baixo
		// Canto inferior esquerdo
		if (params.xLowerBound == 0 && params.yLowerBound == 0)
			idToSend[7] = numberOfProcess - 1;
		// Lateral esquerdo
		else if(params.xLowerBound == 0 && params.yLowerBound != 0)
			idToSend[7] = rank - 1;
		// Baixo
		else if(params.xLowerBound != 0 && params.yLowerBound == 0)
			idToSend[7] = rank + numberOfProcess - params.xSize - 1;
		else
			idToSend[7] = rank - params.xSize - 1;
	}
	// Organizada em linhas
	else {

	}

	// Time Step simulation
	for (k = 0; k < params.timeStep; k++) {
		// Run throw all the cells and resets all the center of mass
		memset(grid.m, 0, params.gridSize*sizeof(double));
		memset(grid.centerOfMassX, 0, params.gridSize*sizeof(double));
		memset(grid.centerOfMassY, 0, params.gridSize*sizeof(double));

		for(int i = params.n_part - 1; i >= 0; i--) {
			if(par[i].number >= 0) {
				CENTEROFMASSX(par[i].gridCoordinateX, par[i].gridCoordinateY) = CENTEROFMASSX(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m * par[i].positionX;
				CENTEROFMASSY(par[i].gridCoordinateX, par[i].gridCoordinateY) = CENTEROFMASSY(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m * par[i].positionY;
				MASS(par[i].gridCoordinateX, par[i].gridCoordinateY) = MASS(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m;
			}
		}
		for (int i = params.yUpperBound; i >= params.yLowerBound; i--) {
			for (int j = params.xUpperBound; j >= params.xLowerBound; j--) {
				CENTEROFMASSX(i, j) = CENTEROFMASSX(i, j)/MASS(i, j);
				CENTEROFMASSY(i, j) = CENTEROFMASSY(i, j)/MASS(i, j);
			}
		}


		// Copies the center of mass to be transmmited
		int pos = 0;
		for (int i = params.yUpperBound; i >= params.yLowerBound; i--) {
			gridSendReceive[LEFTPROCESS][pos].centerOfMassX = CENTEROFMASSX(i, params.xLowerBound);
			gridSendReceive[LEFTPROCESS][pos].centerOfMassY = CENTEROFMASSY(i, params.xLowerBound);
			gridSendReceive[LEFTPROCESS][pos].m = MASS(i, params.xLowerBound);
			gridSendReceive[RIGHTPROCESS][pos].centerOfMassX = CENTEROFMASSX(i, params.xUpperBound);
			gridSendReceive[RIGHTPROCESS][pos].centerOfMassY = CENTEROFMASSY(i, params.xUpperBound);
			gridSendReceive[RIGHTPROCESS][pos++].m = MASS(i, params.xUpperBound);
		}
		
		pos = 0;
		for (int j = params.xUpperBound; j >= params.xLowerBound; j--) {
			gridSendReceive[UPPROCESS][pos].centerOfMassX = CENTEROFMASSX(params.yUpperBound, j);
			gridSendReceive[UPPROCESS][pos].centerOfMassY = CENTEROFMASSY(params.yUpperBound, j);
			gridSendReceive[UPPROCESS][pos].m = MASS(params.yUpperBound, j);
			gridSendReceive[DOWNPROCESS][pos].centerOfMassX = CENTEROFMASSX(params.yLowerBound, j);
			gridSendReceive[DOWNPROCESS][pos].centerOfMassY = CENTEROFMASSY(params.yLowerBound, j);
			gridSendReceive[DOWNPROCESS][pos++].m = MASS(params.yLowerBound, j);
		}
		gridSendReceive[UPLEFTPROCESS][0].centerOfMassX = CENTEROFMASSX(params.yUpperBound, params.xLowerBound);
		gridSendReceive[UPLEFTPROCESS][0].centerOfMassY = CENTEROFMASSY(params.yUpperBound, params.xLowerBound);
		gridSendReceive[UPLEFTPROCESS][0].m = MASS(params.yUpperBound, params.xLowerBound);
		
		gridSendReceive[UPRIGHTPROCESS][0].centerOfMassX = CENTEROFMASSX(params.yUpperBound, params.xUpperBound);
		gridSendReceive[UPRIGHTPROCESS][0].centerOfMassY = CENTEROFMASSY(params.yUpperBound, params.xUpperBound);
		gridSendReceive[UPRIGHTPROCESS][0].m = MASS(params.yUpperBound, params.xUpperBound);

		gridSendReceive[DOWNRIGHTPROCESS][0].centerOfMassX = CENTEROFMASSX(params.yLowerBound, params.xUpperBound);
		gridSendReceive[DOWNRIGHTPROCESS][0].centerOfMassY = CENTEROFMASSY(params.yLowerBound, params.xUpperBound);
		gridSendReceive[DOWNRIGHTPROCESS][0].m = MASS(params.yLowerBound, params.xUpperBound);

		gridSendReceive[DOWNLEFTPROCESS][0].centerOfMassX = CENTEROFMASSX(params.yLowerBound, params.xLowerBound);
		gridSendReceive[DOWNLEFTPROCESS][0].centerOfMassY = CENTEROFMASSY(params.yLowerBound, params.xLowerBound);
		gridSendReceive[DOWNLEFTPROCESS][0].m = MASS(params.yLowerBound, params.xLowerBound);

		// Envia os centros de massa das fronteiras e recebe das regiões vizinhas
		if(params.xSize != 1) {
			MPI_Irecv(gridSendReceive[8], params.sizeBigD, mpi_grid_t, idToSend[0], SENDCENTER,MPI_COMM_WORLD, &request);
			MPI_Irecv(gridSendReceive[9], 1, mpi_grid_t, idToSend[1], SENDCENTER,MPI_COMM_WORLD, &request);
			MPI_Irecv(gridSendReceive[10], params.sizeSmallD, mpi_grid_t, idToSend[2], SENDCENTER,MPI_COMM_WORLD, &request);
			MPI_Irecv(gridSendReceive[11], 1, mpi_grid_t, idToSend[3], SENDCENTER,MPI_COMM_WORLD, &request);
			MPI_Irecv(gridSendReceive[12], params.sizeBigD, mpi_grid_t, idToSend[4], SENDCENTER,MPI_COMM_WORLD, &request);
			MPI_Irecv(gridSendReceive[13], 1, mpi_grid_t, idToSend[5], SENDCENTER,MPI_COMM_WORLD, &request);
			MPI_Irecv(gridSendReceive[14], params.sizeSmallD, mpi_grid_t, idToSend[6], SENDCENTER,MPI_COMM_WORLD, &request);
			MPI_Irecv(gridSendReceive[15], 1, mpi_grid_t, idToSend[7], SENDCENTER,MPI_COMM_WORLD, &request);
			MPI_Send(gridSendReceive[0], params.sizeBigD, mpi_grid_t, idToSend[0], 0, MPI_COMM_WORLD);
			MPI_Send(gridSendReceive[1], 1, mpi_grid_t, idToSend[1], 0, MPI_COMM_WORLD);
			MPI_Send(gridSendReceive[2], params.sizeSmallD, mpi_grid_t, idToSend[2], 0, MPI_COMM_WORLD);
			MPI_Send(gridSendReceive[3], 1, mpi_grid_t, idToSend[3], 0, MPI_COMM_WORLD);
			MPI_Send(gridSendReceive[4], params.sizeBigD, mpi_grid_t, idToSend[4], 0, MPI_COMM_WORLD);
			MPI_Send(gridSendReceive[5], 1, mpi_grid_t, idToSend[5], 0, MPI_COMM_WORLD);
			MPI_Send(gridSendReceive[6], params.sizeSmallD, mpi_grid_t, idToSend[6], 0, MPI_COMM_WORLD);
			MPI_Send(gridSendReceive[7], 1, mpi_grid_t, idToSend[7], 0, MPI_COMM_WORLD);
			MPI_Wait(&request, &status);
		}

		// Updates the new values of the center of mass
		pos = 0;
		for (int i = params.yUpperBound; i >= params.yLowerBound; i--) {
			CENTEROFMASSX(i, params.xLowerBound) = gridSendReceive[LEFTPROCESS + 8][pos].centerOfMassX;
			CENTEROFMASSY(i, params.xLowerBound) = gridSendReceive[LEFTPROCESS + 8][pos].centerOfMassY;
			MASS(i, params.xLowerBound) = gridSendReceive[LEFTPROCESS + 8][pos].m;
			CENTEROFMASSX(i, params.xUpperBound) = gridSendReceive[RIGHTPROCESS + 8][pos].centerOfMassX;
			CENTEROFMASSY(i, params.xUpperBound) = gridSendReceive[RIGHTPROCESS + 8][pos].centerOfMassY;
			MASS(i, params.xUpperBound) = gridSendReceive[RIGHTPROCESS + 8][pos++].m;
		}
		
		pos = 0;
		for (int j = params.xUpperBound; j >= params.xLowerBound; j--) {
			CENTEROFMASSX(params.yUpperBound, j) = gridSendReceive[UPPROCESS + 8][pos].centerOfMassX;
			CENTEROFMASSY(params.yUpperBound, j) = gridSendReceive[UPPROCESS + 8][pos].centerOfMassY;
			MASS(params.yUpperBound, j) = gridSendReceive[UPPROCESS + 8][pos].m;
			CENTEROFMASSX(params.yLowerBound, j) = gridSendReceive[DOWNPROCESS + 8][pos].centerOfMassX;
			CENTEROFMASSY(params.yLowerBound, j) = gridSendReceive[DOWNPROCESS + 8][pos].centerOfMassY;
			MASS(params.yLowerBound, j) = gridSendReceive[DOWNPROCESS + 8][pos++].m;
		}

		CENTEROFMASSX(params.yUpperBound, params.xLowerBound) = gridSendReceive[UPLEFTPROCESS + 8][0].centerOfMassX;
		CENTEROFMASSY(params.yUpperBound, params.xLowerBound) = gridSendReceive[UPLEFTPROCESS + 8][0].centerOfMassY;
		MASS(params.yUpperBound, params.xLowerBound) = gridSendReceive[UPLEFTPROCESS + 8][0].m;
		
		CENTEROFMASSX(params.yUpperBound, params.xUpperBound) = gridSendReceive[UPRIGHTPROCESS + 8][0].centerOfMassX;
		CENTEROFMASSY(params.yUpperBound, params.xUpperBound) = gridSendReceive[UPRIGHTPROCESS + 8][0].centerOfMassY;
		MASS(params.yUpperBound, params.xUpperBound) = gridSendReceive[UPRIGHTPROCESS + 8][0].m;

		CENTEROFMASSX(params.yLowerBound, params.xUpperBound) = gridSendReceive[DOWNRIGHTPROCESS + 8][0].centerOfMassX;
		CENTEROFMASSY(params.yLowerBound, params.xUpperBound) = gridSendReceive[DOWNRIGHTPROCESS + 8][0].centerOfMassY;
		MASS(params.yLowerBound, params.xUpperBound) = gridSendReceive[DOWNRIGHTPROCESS + 8][0].m;

		CENTEROFMASSX(params.yLowerBound, params.xLowerBound) = gridSendReceive[DOWNLEFTPROCESS + 8][0].centerOfMassX;
		CENTEROFMASSY(params.yLowerBound, params.xLowerBound) = gridSendReceive[DOWNLEFTPROCESS + 8][0].centerOfMassY;
		MASS(params.yLowerBound, params.xLowerBound) = gridSendReceive[DOWNLEFTPROCESS + 8][0].m;

		

		// Compute interactions
		// Run all particles
		long aux1, aux2;
		double invM;
		int sideUPDOWN;
		int sideLEFTRIGHT;
		for(int i = params.n_part - 1; i >= 0; i = i - 1){
			if(par[i].number >= 0) {
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

				// Verificar se particula ficou fora da área de trabalho
				if(par[i].gridCoordinateX < params.xLowerBound || par[i].gridCoordinateX > params.xUpperBound || 
					par[i].gridCoordinateY < params.yLowerBound || par[i].gridCoordinateY > params.yUpperBound) {
					long long destiny = par[i].gridCoordinateX + par[i].gridCoordinateY * params.xSize;
					MPI_Isend(&par[i], 1, mpi_particle_t, destiny, 2, MPI_COMM_WORLD, &request);
				}
			}
		}
		MPI_Wait(&request, &status);
		// Verifica se existe alguma mensagem para ser recebida
		MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);
		if(flag == 1) {

		}
	}

	MPI_Type_free(&mpi_particle_t);
	MPI_Finalize();

	return 0;
}

