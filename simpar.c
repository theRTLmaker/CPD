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
#include <omp.h>

#define SENDCENTER 0

/* Fazer um handler para tratar dos argumentos de entrada e passar para o init_particles*/

int main(int argc, char *argv[])
{ 
	int rank, numberOfProcess, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status; 
	MPI_Comm comm; 
	int *idToSend;

	grid_t grid;
	grid_tt **gridSendReceive;
	particle_t *par;
	particle_t_reduced *parReceive;
	particle_t_reduced **parSend;
	long sizeParReceive;
	long sizeParSend;
	long parAuxX;
	long parAuxY;
	int parSendPos[8];
	int count;
	int provided;

	long k;
	int flag, destiny;

	MPI_Request request[8];
	MPI_Status statuss[8];

	MPI_Init_thread( &argc, &argv, 2, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );
	MPI_Get_processor_name(processor_name, &namelen);

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
	const int nitemsPart = 4;
	MPI_Aint displacementsPart[4] = {	offsetof(particle_t_final, number), 
										offsetof(particle_t_final, positionX), 
										offsetof(particle_t_final, positionY),
										offsetof(particle_t_final, m)
									};

	int block_lengthsPart[4]  = {1, 1, 1, 1};
	MPI_Datatype typesPart[4] = {MPI_LONG_LONG_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_FLOAT};
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

	// Find the division of the grid by the processes
	numberOfProcess = findGridDivision(numberOfProcess, rank);

	int isActive = (rank < numberOfProcess);
	MPI_Comm_split(MPI_COMM_WORLD, isActive, rank, &comm);
	
	if(rank < numberOfProcess) {
		par = CreateParticleArray(params.n_part);

		parReceive =  initParReceived(params.n_part, &sizeParReceive);

		parSend = initParSend(params.n_part, &sizeParSend);

		grid = initTotalGrid(grid, params.ncside);

		gridSendReceive = initGridSendReceive(rank);
		
		
		// Inicia as particulas e descobre a sua posicao na grelha (pode ser paralelizado a descoberta )
		init_particles(par);
		
		#pragma omp parallel
		{
			#pragma omp for
			// Ativam as particulas que lhes pertencem
			for(long long i = params.n_part - 1; i >= 0; i = i - 1) {
				if(par[i].gridCoordinateX >= params.xLowerBound && par[i].gridCoordinateX <= params.xUpperBound && 
				   par[i].gridCoordinateY >= params.yLowerBound && par[i].gridCoordinateY <= params.yUpperBound) {
		        	par[i].number = i;
				}
			}	
		}
		if((idToSend = (int *)malloc(8*sizeof(int))) == NULL) {
			printf("ERROR malloc idToSend\n");fflush(stdout);
			exit(0);
		}

		idToSend = findNeighborsRank(idToSend, rank, numberOfProcess);


		printf("rank %d:\n%d - %d - %d\n%d -   - %d\n%d - %d - %d\n", rank, idToSend[1],
            idToSend[2], idToSend[3], idToSend[0], idToSend[4], idToSend[7], idToSend[6], idToSend[5]);fflush(stdout);

		
		
		int pos = 0;
		// Time Step simulation
		for(k = params.timeStep; k > 0; k = k - 1) {

			if(rank == 0) {
				printf("iteration %ld\n", k); fflush(stdout);
			}

			// Clear the memory position to send and receive particles
			memset(parSend[0], 0, sizeParSend*sizeof(particle_t_reduced));
			memset(parSend[1], 0, sizeParSend*sizeof(particle_t_reduced));
			memset(parSend[2], 0, sizeParSend*sizeof(particle_t_reduced));
			memset(parSend[3], 0, sizeParSend*sizeof(particle_t_reduced));
			memset(parSend[4], 0, sizeParSend*sizeof(particle_t_reduced));
			memset(parSend[5], 0, sizeParSend*sizeof(particle_t_reduced));
			memset(parSend[6], 0, sizeParSend*sizeof(particle_t_reduced));
			memset(parSend[7], 0, sizeParSend*sizeof(particle_t_reduced));
			memset(parReceive, 0, sizeParReceive*sizeof(particle_t_reduced));
			memset(parSendPos, 0, 8*sizeof(int));

			// Run throw all the cells and resets all the center of mass
			memset(grid.m, 0, params.gridSize*sizeof(double));
			memset(grid.centerOfMassX, 0, params.gridSize*sizeof(double));
			memset(grid.centerOfMassY, 0, params.gridSize*sizeof(double));

			#pragma omp parallel
			{
				#pragma omp for
				for(int i = params.n_part - 1; i >= 0; i = i - 1) {
					if(par[i].number >= 0) {
						CENTEROFMASSX(par[i].gridCoordinateX, par[i].gridCoordinateY) = CENTEROFMASSX(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m * par[i].positionX;
						CENTEROFMASSY(par[i].gridCoordinateX, par[i].gridCoordinateY) = CENTEROFMASSY(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m * par[i].positionY;
						MASS(par[i].gridCoordinateX, par[i].gridCoordinateY) = MASS(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m;
					}
				}
				#pragma omp for
				for (int i = params.yUpperBound; i >= params.yLowerBound; i = i - 1) {
					for (int j = params.xUpperBound; j >= params.xLowerBound; j = j - 1) {
						CENTEROFMASSX(i, j) = CENTEROFMASSX(i, j)/MASS(i, j);
						CENTEROFMASSY(i, j) = CENTEROFMASSY(i, j)/MASS(i, j);
					}
				}
			}
			
			pos = 0;
			// Copies the center of mass to be transmmited
			for (int i = params.yUpperBound; i >= params.yLowerBound; i = i - 1) {
				gridSendReceive[LEFTPROCESS][pos].centerOfMassX = CENTEROFMASSX(i, params.xLowerBound);
				gridSendReceive[LEFTPROCESS][pos].centerOfMassY = CENTEROFMASSY(i, params.xLowerBound);
				gridSendReceive[LEFTPROCESS][pos].m = MASS(i, params.xLowerBound);
				gridSendReceive[RIGHTPROCESS][pos].centerOfMassX = CENTEROFMASSX(i, params.xUpperBound);
				gridSendReceive[RIGHTPROCESS][pos].centerOfMassY = CENTEROFMASSY(i, params.xUpperBound);
				gridSendReceive[RIGHTPROCESS][pos++].m = MASS(i, params.xUpperBound);
			}

			pos = 0;
			for (int j = params.xUpperBound; j >= params.xLowerBound; j = j- 1) {
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
			MPI_Irecv(gridSendReceive[8], params.sizeVertical, mpi_grid_t, idToSend[0], SENDCENTER,comm, &request[0]);
			MPI_Irecv(gridSendReceive[9], 1, mpi_grid_t, idToSend[1], SENDCENTER,comm, &request[1]);
			MPI_Irecv(gridSendReceive[10], params.sizeHorizontal, mpi_grid_t, idToSend[2], SENDCENTER,comm, &request[2]);
			MPI_Irecv(gridSendReceive[11], 1, mpi_grid_t, idToSend[3], SENDCENTER,comm, &request[3]);
			MPI_Irecv(gridSendReceive[12], params.sizeVertical, mpi_grid_t, idToSend[4], SENDCENTER,comm, &request[4]);
			MPI_Irecv(gridSendReceive[13], 1, mpi_grid_t, idToSend[5], SENDCENTER,comm, &request[5]);
			MPI_Irecv(gridSendReceive[14], params.sizeHorizontal, mpi_grid_t, idToSend[6], SENDCENTER,comm, &request[6]);
			MPI_Irecv(gridSendReceive[15], 1, mpi_grid_t, idToSend[7], SENDCENTER,comm, &request[7]);
		
			MPI_Send(gridSendReceive[0], params.sizeVertical, mpi_grid_t, idToSend[0], 0, comm);
			MPI_Send(gridSendReceive[1], 1, mpi_grid_t, idToSend[1], 0, comm);
			MPI_Send(gridSendReceive[2], params.sizeHorizontal, mpi_grid_t, idToSend[2], 0, comm);
			MPI_Send(gridSendReceive[3], 1, mpi_grid_t, idToSend[3], 0, comm);
			MPI_Send(gridSendReceive[4], params.sizeVertical, mpi_grid_t, idToSend[4], 0, comm);
			MPI_Send(gridSendReceive[6], params.sizeHorizontal, mpi_grid_t, idToSend[6], 0, comm);
			MPI_Send(gridSendReceive[5], 1, mpi_grid_t, idToSend[5], 0, comm);
			MPI_Send(gridSendReceive[7], 1, mpi_grid_t, idToSend[7], 0, comm);

			MPI_Waitall(8, request, statuss);



			// Updates the new values of the center of mass
			pos = 0;
			for (int i = params.yUpperBound; i >= params.yLowerBound; i = i - 1) {
				CENTEROFMASSX(i, params.xLowerBound) = gridSendReceive[LEFTPROCESS + 8][pos].centerOfMassX;
				CENTEROFMASSY(i, params.xLowerBound) = gridSendReceive[LEFTPROCESS + 8][pos].centerOfMassY;
				MASS(i, params.xLowerBound) = gridSendReceive[LEFTPROCESS + 8][pos].m;
				CENTEROFMASSX(i, params.xUpperBound) = gridSendReceive[RIGHTPROCESS + 8][pos].centerOfMassX;
				CENTEROFMASSY(i, params.xUpperBound) = gridSendReceive[RIGHTPROCESS + 8][pos].centerOfMassY;
				MASS(i, params.xUpperBound) = gridSendReceive[RIGHTPROCESS + 8][pos++].m;
			}
			
			pos = 0;
			for (int j = params.xUpperBound; j >= params.xLowerBound; j = j - 1) {
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

			#pragma omp for 
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

					parAuxX = par[i].positionX * params.ncside;
		            parAuxY = par[i].positionY * params.ncside;

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
						if(parAuxX < params.xLowerBound) {
							if(parAuxY < params.yLowerBound)
								destiny = 7;
							else if(parAuxY > params.yUpperBound)
								destiny = 5;
							else
								destiny = 0;
						}
						else if(parAuxX > params.xUpperBound) {
							if(parAuxY < params.yLowerBound)
								destiny = 5;
							else if(parAuxY > params.yUpperBound)
								destiny = 3;
							else
								destiny = 4;
						}
						else{
							if(parAuxY < params.yLowerBound)
								destiny = 6;
							else
								destiny = 2;
						}

						parSend[destiny][parSendPos[destiny]].number = par[i].number;
			            parSend[destiny][parSendPos[destiny]].positionX = par[i].positionX;
			        	parSend[destiny][parSendPos[destiny]].positionY = par[i].positionY;
			            parSend[destiny][parSendPos[destiny]].vx = par[i].vx;
			            parSend[destiny][parSendPos[destiny]].vy = par[i].vy;
			            parSend[destiny][parSendPos[destiny]].gridCoordinateX = par[i].gridCoordinateX;
			            parSend[destiny][parSendPos[destiny]++].gridCoordinateY = par[i].gridCoordinateY;
			            
						par[i].number = -1;
					}
				}
			}
			
			// Send the particles
			for(int i = 0; i < 8; ++i) {
				if(parSendPos[i] != 0) {
					MPI_Isend(parSend[i], parSendPos[i], mpi_particle_t_reduced, idToSend[i], 2 , comm, &request[i]);
				}
			}

			// Barreira de sincronizacao
			if(MPI_Barrier(comm) != MPI_SUCCESS) {
				printf(" Error on barrier on iteration %ld\n", k); fflush(stdout);
			}

 			for (int i = 0; i < 8; ++i) {
	 			do{	
		 			MPI_Iprobe(idToSend[i], 2, comm, &flag, &status);
					if(flag) { 
		            	MPI_Recv(parReceive, sizeParReceive, mpi_particle_t_reduced, idToSend[i], 2, comm, &status);
		            	MPI_Get_count(&status, mpi_particle_t_reduced, &count);
		            	for (int i = 0; i < count; ++i){
		            		par[parReceive[i].number].number = parReceive[i].number;
			            	par[parReceive[i].number].positionX = parReceive[i].positionX;
			            	par[parReceive[i].number].positionY = parReceive[i].positionY;
			            	par[parReceive[i].number].vx = parReceive[i].vx;
			            	par[parReceive[i].number].vy = parReceive[i].vy;
			            	par[parReceive[i].number].gridCoordinateX = parReceive[i].gridCoordinateX;
			            	par[parReceive[i].number].gridCoordinateY = parReceive[i].gridCoordinateY;
		            	}
					}
				}while(flag);
			}

			for (int i = 0; i < 8; ++i) {
				if(parSendPos[i] != 0) {
					MPI_Wait(&request[i], MPI_STATUS_IGNORE);
				}
			} 
		}

		particle_t_final particle_recv;
		if(rank != 0) {
			for(long long i = params.n_part - 1; i >= 0; i = i - 1) {
				if(par[i].number != -1) {
					particle_recv.number = par[i].number;
					particle_recv.positionX = par[i].positionX;
					particle_recv.positionY = par[i].positionY;
					particle_recv.m = par[i].m;
					MPI_Send(&particle_recv, 1, mpi_particle_t, 0, 3, comm);
					if(par[i].number == 0) {
						printf("%.2f %.2f\n", par[i].positionX, par[i].positionY);fflush(stdout);
					}
				}
			}
		}
		else {
			// Computes the total center of mass
			double centerOfMassX = 0;
			double centerOfMassY = 0;
			double totalMass = 0;

			// Calculates the center of mass of all cells
			for(long long i = params.n_part - 1; i >= 0; i = i - 1) {
				if(par[i].number == -1) {
					MPI_Recv(&particle_recv, 1, mpi_particle_t, MPI_ANY_SOURCE, 3, comm, &status);
					par[i].positionX = particle_recv.positionX;
					par[i].positionY = particle_recv.positionY;
					par[i].m = particle_recv.m;
				}
				centerOfMassX = centerOfMassX + par[i].m * par[i].positionX;
				centerOfMassY = centerOfMassY + par[i].m * par[i].positionY;
				totalMass = totalMass + par[i].m;
			}

			// prints the information
			centerOfMassX = centerOfMassX / totalMass;
			centerOfMassY = centerOfMassY / totalMass;
			if(par[0].number == 0) { 
				printf("%.2f %.2f\n", par[0].positionX, par[0].positionY);fflush(stdout);
			}
			printf("%.2f %.2f\n", centerOfMassX, centerOfMassY);
	    }
	}
	// Free everything
	MPI_Type_free(&mpi_grid_t);
	MPI_Type_free(&mpi_particle_t_reduced);
	MPI_Type_free(&mpi_particle_t_complete);
	MPI_Type_free(&mpi_particle_t);
	freeEverything(par, grid, params.ncside);
	free(idToSend);
	MPI_Finalize();

	return 0;
}

