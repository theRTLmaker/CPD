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

	// Time variables
	double startIter = 0;
	double endIter = 0;
	double timeRearrange = 0;
	double timeSending = 0;
	double timeBarrier = 0;
	double timeReciving = 0;
	double start = 0;
	double end = 0;

	grid_t grid;
	grid_tt **gridSendReceive;
	particle_t *par;
	particle_t_reduced *parReceive;
	particle_t_reduced **parSend;
	long sizeParReceive = 0;
	long sizeParSend[8] = {0};
	long incSizeParReceive = 0;
	long incSizeParSend = 0;
	long parAuxX;
	long parAuxY;
	int parSendPos[8];
	int count;
	int provided;

	long k;
	int flag, destiny;

	MPI_Request request[8];
	MPI_Status statuss[8];

	MPI_Init_thread( &argc, &argv, 1, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcess);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );
	MPI_Get_processor_name(processor_name, &namelen);

	// Criação de estrutura particula para enviar em MPI
	const int nitemsPart = 3;
	MPI_Aint displacementsPart[3] = {	offsetof(particle_t_final, positionX), 
										offsetof(particle_t_final, positionY),
										offsetof(particle_t_final, m)
									};

	int block_lengthsPart[3]  = {1, 1, 1};
	MPI_Datatype typesPart[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype mpi_particle_t_final;
	MPI_Type_create_struct(nitemsPart, block_lengthsPart, displacementsPart, 
							typesPart, &mpi_particle_t_final);
	MPI_Type_commit(&mpi_particle_t_final);

	// Criação de estrutura particula para enviar em MPI
	const int nitemsPartReduce = 8;
	MPI_Aint displacementsPartReduced[8] = {	offsetof(particle_t_reduced, isZero), 
												offsetof(particle_t_reduced, m), 
												offsetof(particle_t_reduced, positionX), 
												offsetof(particle_t_reduced, positionY),
												offsetof(particle_t_reduced, vx),
												offsetof(particle_t_reduced, vy),
												offsetof(particle_t_reduced, gridCoordinateX),
												offsetof(particle_t_reduced, gridCoordinateY)
											};

	int block_lengthsPartReduced[8]  = {1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype typesPartReduced[8] = {MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_LONG};
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

	start = MPI_Wtime();

	// init particles
	handler_input(argc, argv);

	// Find the division of the grid by the processes
	numberOfProcess = findGridDivision(numberOfProcess, rank);

	int isActive = (rank < numberOfProcess);
	MPI_Comm_split(MPI_COMM_WORLD, isActive, rank, &comm);
	
	if(rank < numberOfProcess) {
		// Cria um array de particulas com dimensao 2*(n_part/numberOfProcess)
		par = CreateParticleArray(numberOfProcess);
		// Inicia as particulas que estao na zona da grelha de controlo
		par = init_particles(par, numberOfProcess, rank);
		//printf("rank %d - partVectSize %lld, activeParticles %lld\n", rank, params.partVectSize, params.activeParticles);fflush(stdout);
		parReceive =  initParReceived(params.n_part, &sizeParReceive, rank, &incSizeParReceive);
		/*if(rank == 0) {
			printf("rank %d - sizeParReceive %ld, incSizeParReceive %ld\n", rank, sizeParReceive, incSizeParReceive);fflush(stdout);
		}*/
		
		parSend = initParSend(params.n_part, sizeParSend, rank, &incSizeParSend);
		/*if(rank == 0) {
			printf("rank %d - sizeParSend %ld, incSizeParSend %ld\n", rank, sizeParSend[0], incSizeParSend);fflush(stdout);
		}*/

		grid = initTotalGrid(grid, params.ncside);

		gridSendReceive = initGridSendReceive(rank);
		
		if((idToSend = (int *)malloc(8*sizeof(int))) == NULL) {
			printf("ERROR malloc idToSend\n");fflush(stdout);
			exit(0);
		}

		idToSend = findNeighborsRank(idToSend, rank, numberOfProcess);


		/*printf("rank %d - %lld:\n%d - %d - %d\n%d -   - %d\n%d - %d - %d\n", rank, params.activeParticles, idToSend[1],
			idToSend[2], idToSend[3], idToSend[0], idToSend[4], idToSend[7], idToSend[6], idToSend[5]);fflush(stdout);*/

		startIter = MPI_Wtime();
		
		int pos = 0;
		// Time Step simulation
		for(k = params.timeStep; k > 0; k = k - 1) {

			if(rank == 0) {
				printf("iteration %ld\n", k); fflush(stdout);
			}

			// Clear the memory position to send and receive particles
			memset(parSend[0], 0, sizeParSend[0]*sizeof(particle_t_reduced));
			memset(parSend[1], 0, sizeParSend[1]*sizeof(particle_t_reduced));
			memset(parSend[2], 0, sizeParSend[2]*sizeof(particle_t_reduced));
			memset(parSend[3], 0, sizeParSend[3]*sizeof(particle_t_reduced));
			memset(parSend[4], 0, sizeParSend[4]*sizeof(particle_t_reduced));
			memset(parSend[5], 0, sizeParSend[5]*sizeof(particle_t_reduced));
			memset(parSend[6], 0, sizeParSend[6]*sizeof(particle_t_reduced));
			memset(parSend[7], 0, sizeParSend[7]*sizeof(particle_t_reduced));
			memset(parReceive, 0, sizeParReceive*sizeof(particle_t_reduced));
			memset(parSendPos, 0, 8*sizeof(int));

			// Run throw all the cells and resets all the center of mass
			memset(grid.m, 0, params.gridSize*sizeof(double));
			memset(grid.centerOfMassX, 0, params.gridSize*sizeof(double));
			memset(grid.centerOfMassY, 0, params.gridSize*sizeof(double));

			#pragma omp parallel
			{
				#pragma omp for
				for(int i = params.activeParticles - 1; i >= 0; i = i - 1) {
					CENTEROFMASSX(par[i].gridCoordinateX, par[i].gridCoordinateY) = CENTEROFMASSX(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m * par[i].positionX;
					CENTEROFMASSY(par[i].gridCoordinateX, par[i].gridCoordinateY) = CENTEROFMASSY(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m * par[i].positionY;
					MASS(par[i].gridCoordinateX, par[i].gridCoordinateY) = MASS(par[i].gridCoordinateX, par[i].gridCoordinateY) + par[i].m;
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
				gridSendReceive[RIGHTPROCESS][pos].m = MASS(i, params.xUpperBound);
				pos = pos + 1;
			}

			pos = 0;
			for (int j = params.xUpperBound; j >= params.xLowerBound; j = j- 1) {
				gridSendReceive[UPPROCESS][pos].centerOfMassX = CENTEROFMASSX(params.yUpperBound, j);
				gridSendReceive[UPPROCESS][pos].centerOfMassY = CENTEROFMASSY(params.yUpperBound, j);
				gridSendReceive[UPPROCESS][pos].m = MASS(params.yUpperBound, j);
				gridSendReceive[DOWNPROCESS][pos].centerOfMassX = CENTEROFMASSX(params.yLowerBound, j);
				gridSendReceive[DOWNPROCESS][pos].centerOfMassY = CENTEROFMASSY(params.yLowerBound, j);
				gridSendReceive[DOWNPROCESS][pos].m = MASS(params.yLowerBound, j);
				pos = pos + 1;
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
				MASS(i, params.xUpperBound) = gridSendReceive[RIGHTPROCESS + 8][pos].m;
				pos = pos + 1;
			}
			
			pos = 0;
			for (int j = params.xUpperBound; j >= params.xLowerBound; j = j - 1) {
				CENTEROFMASSX(params.yUpperBound, j) = gridSendReceive[UPPROCESS + 8][pos].centerOfMassX;
				CENTEROFMASSY(params.yUpperBound, j) = gridSendReceive[UPPROCESS + 8][pos].centerOfMassY;
				MASS(params.yUpperBound, j) = gridSendReceive[UPPROCESS + 8][pos].m;
				CENTEROFMASSX(params.yLowerBound, j) = gridSendReceive[DOWNPROCESS + 8][pos].centerOfMassX;
				CENTEROFMASSY(params.yLowerBound, j) = gridSendReceive[DOWNPROCESS + 8][pos].centerOfMassY;
				MASS(params.yLowerBound, j) = gridSendReceive[DOWNPROCESS + 8][pos].m;
				pos = pos + 1;
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

			timeSending -= MPI_Wtime();

			#pragma omp for 
			for(long long i = params.activeParticles - 1; i >= 0; i = i - 1){
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
				par[i].vx = par[i].vx + par[i].appliedForceX * invM; //a = F/m 
				par[i].vy = par[i].vy + par[i].appliedForceY * invM;

				par[i].positionX = par[i].positionX + par[i].vx + 0.5 * par[i].appliedForceX * invM;//x = x0 + v0t + 0.5 a t^2 (t = 1)
				par[i].positionY = par[i].positionY + par[i].vy + 0.5 * par[i].appliedForceY * invM;

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

					parSend[destiny][parSendPos[destiny]].isZero = par[i].isZero;
					parSend[destiny][parSendPos[destiny]].m = par[i].m;
					parSend[destiny][parSendPos[destiny]].positionX = par[i].positionX;
					parSend[destiny][parSendPos[destiny]].positionY = par[i].positionY;
					parSend[destiny][parSendPos[destiny]].vx = par[i].vx;
					parSend[destiny][parSendPos[destiny]].vy = par[i].vy;
					parSend[destiny][parSendPos[destiny]].gridCoordinateX = par[i].gridCoordinateX;
					parSend[destiny][parSendPos[destiny]].gridCoordinateY = par[i].gridCoordinateY;

					parSendPos[destiny] = parSendPos[destiny] + 1;

					// Caso esgote o espaço, incrementa o tamanho desse vetor
					if(parSendPos[destiny] + 1 >= sizeParSend[destiny]) {
						//printf("--rank %d - parSendPos[%d] %ld incSizeParReceive %ld\n", rank, destiny, sizeParSend[destiny], incSizeParSend);fflush(stdout);
						sizeParSend[destiny] = sizeParSend[destiny] + incSizeParSend;
						if((parSend[destiny] = (particle_t_reduced *)realloc(parSend[destiny], sizeParSend[destiny]*sizeof(particle_t_reduced))) == NULL) {
							printf("ERROR realloc parSend\n");fflush(stdout);
							exit(0);
						}
					}
					
					par[i].active = 0;
				}
			}

			timeSending += MPI_Wtime();

			// reagrupa as particulas todas, retirando o espaço deixado pelas que sairam
			long long underEvaluation = params.activeParticles - 1;
			for (long long i = 0; i < underEvaluation; ++i) {
				timeRearrange -= MPI_Wtime();
				if(par[i].active == 0) {
					while(par[underEvaluation].active == 0 && underEvaluation > i) {
						underEvaluation--;
					}
					if(par[underEvaluation].active != 0) {
						par[i].isZero = par[underEvaluation].isZero;
						par[i].active = par[underEvaluation].active;
						par[i].m = par[underEvaluation].m;
						par[i].positionX = par[underEvaluation].positionX;
						par[i].positionY = par[underEvaluation].positionY;
						par[i].vx = par[underEvaluation].vx;
						par[i].vy = par[underEvaluation].vy;
						par[i].gridCoordinateX = par[underEvaluation].gridCoordinateX;
						par[i].gridCoordinateY = par[underEvaluation].gridCoordinateY;
						par[i].appliedForceX = par[underEvaluation].appliedForceX;
						par[i].appliedForceY = par[underEvaluation].appliedForceY;
						par[underEvaluation].active = 0;
					}
				}
				timeRearrange += MPI_Wtime();
			}

			
			// Send the particles to the adjacent processes
			for(int i = 0; i < 8; ++i) {
				if(parSendPos[i] != 0) {
					timeSending -= MPI_Wtime();
					//printf("\trank %d sent %d particle to %d\n", rank, parSendPos[i], idToSend[i]); fflush(stdout);
					MPI_Isend(parSend[i], parSendPos[i], mpi_particle_t_reduced, idToSend[i], 2 , comm, &request[i]);
					params.activeParticles = params.activeParticles - parSendPos[i];
					timeSending += MPI_Wtime();
				}
			}

			timeBarrier -= MPI_Wtime();
			// Barreira de sincronizacao
			if(MPI_Barrier(comm) != MPI_SUCCESS) {
				printf(" Error on barrier on iteration %ld\n", k); fflush(stdout);
			}
			timeBarrier += MPI_Wtime();

			// Recebe as particulas dos outros processos
			for (int i = 0; i < 8; ++i) {
				timeReciving -= MPI_Wtime();
				do{
					MPI_Iprobe(idToSend[i], 2, comm, &flag, &status);
					if(flag) { 
						MPI_Get_count(&status, mpi_particle_t_reduced, &count);

						// Verifica se tem espaço suficiente para receber todas as particulas
						if(count > sizeParReceive) {
							// Incrementa o tamanho até ser suficiente para receber tudo
							while(count > sizeParReceive) {
								sizeParReceive = sizeParReceive + incSizeParReceive;
								//printf("-rank %d - count %d, sizeParReceive %ld incSizeParReceive %ld\n", rank, count, sizeParReceive, incSizeParReceive);fflush(stdout);
							}
							if((parReceive = (particle_t_reduced *)realloc(parReceive, sizeParReceive*sizeof(particle_t_reduced))) == NULL) {
								printf("ERROR realloc\n");fflush(stdout);
								exit(0);
							}
						}

						MPI_Recv(parReceive, sizeParReceive, mpi_particle_t_reduced, idToSend[i], 2, comm, &status);

						// Verifica se tem espaço para guardar as novas particulas
						if(params.activeParticles + count > params.partVectSize) {
							// Caso se esgote o tamanho, aloca mais uma parcela de numero de particulas/processos
							params.partVectSize = params.partVectSize + params.reallocInc;
							// Realoca vetor com espaço necessario
							if((par = (particle_t *)realloc(par, params.partVectSize*sizeof(particle_t))) == NULL) {
								printf("ERROR malloc\n");fflush(stdout);
								exit(0);
							}
						}
						// Coloca as novas particulas no sitio correto
						for (int j = 0; j < count; ++j){
							par[params.activeParticles].isZero 		= parReceive[j].isZero;
							par[params.activeParticles].m 			= parReceive[j].m;
							par[params.activeParticles].positionX 	= parReceive[j].positionX;
							par[params.activeParticles].positionY	= parReceive[j].positionY;
							par[params.activeParticles].vx 			= parReceive[j].vx;
							par[params.activeParticles].vy 			= parReceive[j].vy;
							par[params.activeParticles].gridCoordinateX = parReceive[j].gridCoordinateX;
							par[params.activeParticles].gridCoordinateY = parReceive[j].gridCoordinateY;
							par[params.activeParticles].active = 1;

							params.activeParticles = params.activeParticles + 1;
						}
					}
				}while(flag);
				timeReciving += MPI_Wtime();
			}
			for (int i = 0; i < 8; ++i) {
				if(parSendPos[i] != 0) {
					MPI_Wait(&request[i], MPI_STATUS_IGNORE);
				}
			}
		}

		endIter = MPI_Wtime();

		// Liberta espaço para novas alocações
		freeParReceive(parReceive);
		freeParSend(parSend);
		freeGridSendReceive(gridSendReceive);
		freeGrid(grid);
		free(idToSend);

		// Computes the total center of mass of its process
		double centerOfMassX = 0;
		double centerOfMassY = 0;
		double totalMass = 0;

		grid_tt final_send;
		final_send.centerOfMassX = 0;
		final_send.centerOfMassY = 0;
		final_send.m = 0;
		if(rank != 0) {
			// Calcula o centro de massa com as particulas do processo 0
			for(long long i = params.activeParticles - 1; i >= 0; i = i - 1) {	
				final_send.centerOfMassX = final_send.centerOfMassX + par[i].m * par[i].positionX;
				final_send.centerOfMassY = final_send.centerOfMassY + par[i].m * par[i].positionY;
				final_send.m = final_send.m + par[i].m;
				if(par[i].isZero != 0) {
					printf("%.2f %.2f\n", par[i].positionX, par[i].positionY);fflush(stdout);
				}
			}
			freeParticles(par);

			// Envia particulas para o processo 0
			MPI_Send(&final_send, 1, mpi_grid_t, 0, 3, comm);
		}
		else {
			// Calcula o centro de massa com as particulas do processo 0
			for(long long i = params.activeParticles - 1; i >= 0; i = i - 1) {				
				centerOfMassX = centerOfMassX + par[i].m * par[i].positionX;
				centerOfMassY = centerOfMassY + par[i].m * par[i].positionY;
				totalMass = totalMass + par[i].m;
				if(par[i].isZero != 0) {
					printf("%.2f %.2f\n", par[i].positionX, par[i].positionY);fflush(stdout);
				}
			}

			// Liberto as particulas para ter espaço para alocar mais
			freeParticles(par);

			// Receve particulas dos demais processos
			for (int i = 1; i < numberOfProcess; ++i) {
				MPI_Recv(&final_send, 1, mpi_grid_t, i, 3, comm, &status);
				centerOfMassX = centerOfMassX + final_send.centerOfMassX;
				centerOfMassY = centerOfMassY + final_send.centerOfMassY;
				totalMass = totalMass + final_send.m;
			}

			// Imprime a informacao do centro de massa
			centerOfMassX = centerOfMassX / totalMass;
			centerOfMassY = centerOfMassY / totalMass;
			printf("%.2f %.2f\n", centerOfMassX, centerOfMassY);fflush(stdout);
		}


		// Barreira de sincronizacao
		if(MPI_Barrier(comm) != MPI_SUCCESS) {
			printf(" Error on barrier on iteration %ld\n", k); fflush(stdout);
		}

		end = MPI_Wtime();
	}

	printf("rank %d - total %.2f seg, Iter %.2f seg, Arrange %.2f seg\n\t\t\ttimeSending %.2f seg, timeReciving %.2f seg\n\t\t\ttimeBarrier %.2f seg\n", rank, end - start, endIter - startIter, timeRearrange, timeSending, timeReciving, timeBarrier);
	// Free everything
	MPI_Type_free(&mpi_grid_t);
	MPI_Type_free(&mpi_particle_t_reduced);
	MPI_Type_free(&mpi_particle_t_final);
	MPI_Finalize();

	return 0;
}

