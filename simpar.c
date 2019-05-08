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
	MPI_Status status;

	grid_t grid;
	grid_tt **gridSendReceive;
	particle_t *par;
	particle_t_reduced par_aux;
	particle_t_reduced par_aux_send;

	long k;
	int idToSend[8] = {0};
	int flag;
	int count;

	MPI_Request request;

	MPI_Init( &argc, &argv );
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcess);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// Criação de estrutura particula para enviar em MPI
	const int nitemsPartcomplete = 10;
	MPI_Aint displacementsPartcomplete[10] = {	offsetof(particle_t, number), 
												offsetof(particle_t, m),
												offsetof(particle_t, positionX), 
												offsetof(particle_t, positionY),
												offsetof(particle_t, vx),
												offsetof(particle_t, vy),
												offsetof(particle_t, gridCoordinateX),
												offsetof(particle_t, gridCoordinateY),
												offsetof(particle_t, appliedForceX),
												offsetof(particle_t, appliedForceY)
											};

	int block_lengthsPartcomplete[10]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype typesPartcomplete[10] = {	MPI_LONG_LONG_INT, 
											MPI_FLOAT, 
											MPI_DOUBLE, 
											MPI_DOUBLE, 
											MPI_DOUBLE, 
											MPI_DOUBLE, 
											MPI_LONG, 
											MPI_LONG, 
											MPI_DOUBLE, 
											MPI_DOUBLE
										};
	MPI_Datatype mpi_particle_t_complete;
	MPI_Type_create_struct(nitemsPartcomplete, block_lengthsPartcomplete, displacementsPartcomplete, 
							typesPartcomplete, &mpi_particle_t_complete);
	MPI_Type_commit(&mpi_particle_t_complete);

	// Criação de estrutura particula para enviar em MPI
	const int nitemsPart = 3;
	MPI_Aint displacementsPart[3] = {	offsetof(particle_t_reduced, number), 
										offsetof(particle_t_reduced, positionX), 
										offsetof(particle_t_reduced, positionY)
									};

	int block_lengthsPart[3]  = {1, 1, 1};
	MPI_Datatype typesPart[3] = {MPI_LONG_LONG_INT, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype mpi_particle_t;
	MPI_Type_create_struct(nitemsPart, block_lengthsPart, displacementsPart, 
							typesPart, &mpi_particle_t);
	MPI_Type_commit(&mpi_particle_t);

	// Criação de estrutura particula para enviar em MPI
	const int nitemsPartReduce = 7;
	MPI_Aint displacementsPartReduced[7] = {	offsetof(particle_t_reduced, number), 
												offsetof(particle_t_reduced, positionX), 
												offsetof(particle_t_reduced, positionY),
												offsetof(particle_t_reduced, vx),
												offsetof(particle_t_reduced, vy),
												offsetof(particle_t_reduced, gridCoordinateX),
												offsetof(particle_t_reduced, gridCoordinateY)
											};

	int block_lengthsPartReduced[7]  = {1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype typesPartReduced[7] = {MPI_LONG_LONG_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_LONG};
	MPI_Datatype mpi_particle_t_reduced;
	MPI_Type_create_struct(nitemsPartReduce, block_lengthsPartReduced, displacementsPartReduced, 
							typesPartReduced, &mpi_particle_t_reduced);
	MPI_Type_commit(&mpi_particle_t_reduced);

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

	numberOfProcess = findGridDivision(numberOfProcess, rank);
	
	grid = initTotalGrid(grid, params.ncside);

	gridSendReceive = initGridSendReceive(rank);
	
	if(rank == 0) {
		// Inicia as particulas e descobre a sua posicao na grelha (pode ser paralelizado a descoberta )
		init_particles(par);
	}
	
	MPI_Bcast(par, params.n_part, mpi_particle_t_complete, 0,  MPI_COMM_WORLD);

	// Ativam as particulas que lhes pertencem
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
	/*  						topo
	lateral esquerdo	-------------  lateral direito
				    	| 6 | 7 | 8 |
					    -------------
						| 3 | 4 | 5 |
						-------------
						| 0 | 1 | 2 |
						-------------
							baixo

*/


	/*  						topo
	lateral esquerdo	-------------  lateral direito
				    2	|   |       |
					    - 2 -   3   -
					1	|   |       |
						-------------
					0	| 0 |   1   |
						-------------
						  0   1   2
							baixo

*/

	// Organizado em grelha
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

	/* printf("rank %d:\n%d - %d - %d\n%d -   - %d\n%d - %d - %d\n", rank, idToSend[1],
		idToSend[2], idToSend[3], idToSend[0], idToSend[4], idToSend[7], idToSend[6], idToSend[5]);fflush(stdout);
	*/

	// Barreira de sincronizacao
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf( "Error on barrier on iteration\n");fflush(stdout);
	}

	
	int pos = 0;
	// Time Step simulation
	for (k = 0; k < params.timeStep; k++) {

		printf("iteration %ld\n", k);fflush(stdout);


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

		pos = 0;
		// Copies the center of mass to be transmmited
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
		MPI_Irecv(gridSendReceive[8], params.sizeVertical, mpi_grid_t, idToSend[0], SENDCENTER,MPI_COMM_WORLD, &request);
		MPI_Irecv(gridSendReceive[9], 1, mpi_grid_t, idToSend[1], SENDCENTER,MPI_COMM_WORLD, &request);
		//printf("1 - olaaaa %d\n", rank);fflush(stdout);
		MPI_Irecv(gridSendReceive[10], params.sizeHorizontal, mpi_grid_t, idToSend[2], SENDCENTER,MPI_COMM_WORLD, &request);
		//printf("2 - olaaaa %d\n", rank);fflush(stdout);
		MPI_Irecv(gridSendReceive[11], 1, mpi_grid_t, idToSend[3], SENDCENTER,MPI_COMM_WORLD, &request);
		//printf("3 - olaaaa %d\n", rank);fflush(stdout);
		MPI_Irecv(gridSendReceive[12], params.sizeVertical, mpi_grid_t, idToSend[4], SENDCENTER,MPI_COMM_WORLD, &request);
		MPI_Irecv(gridSendReceive[13], 1, mpi_grid_t, idToSend[5], SENDCENTER,MPI_COMM_WORLD, &request);
		MPI_Irecv(gridSendReceive[14], params.sizeHorizontal, mpi_grid_t, idToSend[6], SENDCENTER,MPI_COMM_WORLD, &request);
		MPI_Irecv(gridSendReceive[15], 1, mpi_grid_t, idToSend[7], SENDCENTER,MPI_COMM_WORLD, &request);
		MPI_Send(gridSendReceive[0], params.sizeVertical, mpi_grid_t, idToSend[0], 0, MPI_COMM_WORLD);
		MPI_Send(gridSendReceive[1], 1, mpi_grid_t, idToSend[1], 0, MPI_COMM_WORLD);
		MPI_Send(gridSendReceive[2], params.sizeHorizontal, mpi_grid_t, idToSend[2], 0, MPI_COMM_WORLD);
		MPI_Send(gridSendReceive[3], 1, mpi_grid_t, idToSend[3], 0, MPI_COMM_WORLD);
		MPI_Send(gridSendReceive[4], params.sizeVertical, mpi_grid_t, idToSend[4], 0, MPI_COMM_WORLD);
		MPI_Send(gridSendReceive[5], 1, mpi_grid_t, idToSend[5], 0, MPI_COMM_WORLD);
		MPI_Send(gridSendReceive[6], params.sizeHorizontal, mpi_grid_t, idToSend[6], 0, MPI_COMM_WORLD);
		MPI_Send(gridSendReceive[7], 1, mpi_grid_t, idToSend[7], 0, MPI_COMM_WORLD);
		MPI_Wait(&request, &status);


		// Barreira de sincronizacao
		if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
			printf( "Error on barrier on iteration %ld\n", k);fflush(stdout);
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
long long sent = 0;

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
					par_aux_send.number = par[i].number;
		            par_aux_send.positionX = par[i].positionX;
		            par_aux_send.positionY = par[i].positionY;
		            par_aux_send.vx = par[i].vx;
		            par_aux_send.vy = par[i].vy;
		            par_aux_send.gridCoordinateX = par[i].gridCoordinateX;
		            par_aux_send.gridCoordinateY = par[i].gridCoordinateY;

					int destiny = (par[i].gridCoordinateX / (params.ncside / params.xSize));
					if(destiny == params.xSize) 
						destiny -= 1;
					int destiny_aux = par[i].gridCoordinateY / (params.ncside / params.ySize);
					if(destiny_aux == params.ySize)
						destiny_aux -= 1;

					destiny = destiny + (destiny_aux * params.xSize);
					MPI_Isend(&par_aux_send, 1, mpi_particle_t_reduced, destiny, 2, MPI_COMM_WORLD, &request);
					par[i].number = -1;
					sent++;

					printf("\trank %d sent particle %d , %d-> x: %ld, y: %ld\n", rank, destiny, par[i].gridCoordinateY % params.ySize, par[i].gridCoordinateX, par[i].gridCoordinateY);fflush(stdout);
					

				}
				MPI_Wait(&request, &status);
			}
		}

		// Barreira de sincronizacao
		if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
			printf( "Error on barrier on iteration %ld\n", k);fflush(stdout);
		}

		//printf("rank %d After barrier, iteration %ld\n", rank, k);fflush(stdout);

		// Verifica se existe alguma mensagem para ser recebida
		MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);
long long reeceived = 0;
		//printf("rank %d flag %d, barrier, iteration %ld\n", rank, flag, k);fflush(stdout);
		if(flag == 1) {
			MPI_Get_count( &status, mpi_particle_t_reduced, &count );
			for (long i = 0; i < 8; ++i) {
				MPI_Iprobe(idToSend[i], 2, MPI_COMM_WORLD, &flag, &status);
				if(flag == 1) {
					MPI_Get_count( &status, mpi_particle_t_reduced, &count );
		            if (count != MPI_UNDEFINED) {
			            for (int j = 0; j < count; ++j) {
			            	MPI_Recv(&par_aux, 1, mpi_particle_t_reduced, idToSend[i], 2, MPI_COMM_WORLD, &status);
			            	par[par_aux.number].number = par_aux.number;
			            	par[par_aux.number].positionX = par_aux.positionX;
			            	par[par_aux.number].positionY = par_aux.positionY;
			            	par[par_aux.number].vx = par_aux.vx;
			            	par[par_aux.number].vy = par_aux.vy;
			            	par[par_aux.number].gridCoordinateX = par_aux.gridCoordinateX;
			            	par[par_aux.number].gridCoordinateY = par_aux.gridCoordinateY;
			            	printf("\trank %d ŕeceived particle %d\n", rank, idToSend[i]);fflush(stdout);
			            	reeceived++;
			            }
		            }
				}
			}
		}
		printf("\t\trank %d end iteration %ld -> received %lld, sent %lld\n", rank, k, reeceived, sent);fflush(stdout);
reeceived = 0;			
		for(long long i = 0; i < params.n_part; i++) {
			if(par[i].number != -1) {
				reeceived++;
			}
		}


		printf("\t\trank %d end iteration %ld -> particels %lld\n", rank, k, reeceived);fflush(stdout);
	}


//	printf("rank %d ENDDDDDDDDDDD\n", rank);fflush(stdout);
	// Barreira de sincronizacao
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf( "Error on barrier on iteration %ld\n", k);fflush(stdout);
	}


	printf("\n");fflush(stdout);

	particle_t_final particle_recv;
long long max = 0;
long long reeceived = 0;
	if(rank != 0) {
		for(long long i = 0; i < params.n_part; i++) {
			//printf("rank: %d, i = %lld\n",rank, i);fflush(stdout);
			if(par[i].number != -1) {
				particle_recv.number = par[i].number;
				particle_recv.positionX = par[i].positionX;
				particle_recv.positionY = par[i].positionY;
				MPI_Isend(&particle_recv, 1, mpi_particle_t, 0, 3, MPI_COMM_WORLD, &request);
				max++;
			}
		}
		printf("\t\t\trank: %d -> Send ALL -> max %lld\n", rank, max);fflush(stdout);
	}
	else {
		for(long long i = 0; i < params.n_part; i++) {
			if(par[i].number == -1) {
				//printf("i = %lld -> reeceived %lld\n", i, reeceived);fflush(stdout);
				MPI_Recv(&particle_recv, 1, mpi_particle_t, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
				par[particle_recv.number].positionX = particle_recv.positionX;
				par[particle_recv.number].positionY = particle_recv.positionY;
				reeceived++;
			}
		}

		printf("\t\t\tReceived ALL\n");fflush(stdout);

		// Computes the total center of mass
		double centerOfMassX = 0;
		double centerOfMassY = 0;
		double totalMass = 0;

		// Calculates the center of mass of all cells
		for(long long i = 0; i < params.n_part; i++) {
			centerOfMassX = centerOfMassX + par[i].m * par[i].positionX;
			centerOfMassY = centerOfMassY + par[i].m * par[i].positionY;
			totalMass = totalMass + par[i].m;
		}

		// prints the information
		centerOfMassX = centerOfMassX / totalMass;
		centerOfMassY = centerOfMassY / totalMass;
		printf("%.2f %.2f\n", par[0].positionX, par[0].positionY);
		printf("%.2f %.2f\n", centerOfMassX, centerOfMassY);
    }

	// Free everything
	MPI_Type_free(&mpi_grid_t);
	MPI_Type_free(&mpi_particle_t_reduced);
	MPI_Type_free(&mpi_particle_t_complete);
	MPI_Type_free(&mpi_particle_t);
	freeEverything(par, grid, params.ncside);

	MPI_Finalize();

	return 0;
}

