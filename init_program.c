#include "init_program.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define RND0_1 ((double) random() / ((long long)1<<31))

parameters params; 

void handler_input(int argc ,char *argv[]) {
	if(argc < 5){
		printf("Wrong number of input parameteres\n");
		exit(0);
	}

	params.seed = atol(argv[1]);
	//printf("seed = %ld\n", seed );
	params.ncside = atol(argv[2]);
	//printf("ncside = %ld\n", ncside );
	params.n_part = atoll(argv[3]);
	//printf("n_part = %lld\n", n_part);
	params.timeStep = atoll(argv[4]);

	params.gridSize = params.ncside*params.ncside;

	if(params.seed < 0 || params.ncside < 0 || params.n_part < 0){
		printf("Wrong parameteres\n");
		exit(0);
	}
}

particle_t * CreateParticleArray(int numberOfProcess) {
	params.partVectSize = 2*(params.n_part/numberOfProcess);
	particle_t *par = (particle_t *) malloc(params.partVectSize*sizeof(particle_t));
	if(par == NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	return par;
}

particle_t * init_particles(particle_t *par, int numberOfProcess, int rank) {
    long long i;
    
    srandom(params.seed);

    params.activeParticles = 0;
    
    for(i = 0; i < params.n_part; i++) {    

		par[params.activeParticles].positionX = RND0_1;
		par[params.activeParticles].positionY = RND0_1;

		par[params.activeParticles].vx = RND0_1 / params.ncside / 10.0;
		par[params.activeParticles].vy = RND0_1 / params.ncside / 10.0;

		par[params.activeParticles].m = RND0_1 * params.ncside / (G * 1e6 * params.n_part);

		par[params.activeParticles].gridCoordinateX = par[params.activeParticles].positionX * params.ncside / 1;
		par[params.activeParticles].gridCoordinateY = par[params.activeParticles].positionY * params.ncside / 1;

        // Verifica se particula lhe pertence
        if(par[params.activeParticles].gridCoordinateX >= params.xLowerBound && par[params.activeParticles].gridCoordinateX <= params.xUpperBound && 
        	par[params.activeParticles].gridCoordinateY >= params.yLowerBound && par[params.activeParticles].gridCoordinateY <= params.yUpperBound) {
    		//printf("\t%lld - rank %d x=%ld, y=%ld\n", i, rank, gX, gY);fflush(stdout);
			par[params.activeParticles].active = 1;

			par[params.activeParticles].isZero = 0;
        	
	    	if(i == 0) par[i].isZero = 1;

        	// Incrementa o número de particulas ativas
        	params.activeParticles = params.activeParticles + 1;
    		// Caso se esgote o tamanho, aloca mais uma parcela de numero de particulas/processos
        	if(params.activeParticles + 1 >= params.partVectSize) {
        		params.partVectSize = params.partVectSize + (params.n_part/numberOfProcess);
        		if((par = (particle_t *)realloc((void *) par, params.partVectSize)) == NULL) {
					printf("ERROR malloc\n");
					exit(0);
				}
        	}

        }
    }

    return(par);
}

int findGridDivision(int numberOfProcess, int rank) {
	int dims[2] = {0,0};

	// Se existir mais processes do que celulas, reduz-se
	if(numberOfProcess > params.ncside*params.ncside)
		numberOfProcess = params.ncside*params.ncside;
	// Calcular a divisao da grelha pelos process
	MPI_Dims_create(numberOfProcess, 2, dims);
	if(dims[1] == 1) {
		do {
			numberOfProcess--;
			MPI_Dims_create(numberOfProcess, 2, dims);
		}while(dims[1] == 1);
	}
	params.xSize = dims[1];
	params.ySize = dims[0];

	params.sizeHorizontal = params.ncside / params.xSize;
	long sizeHorizontalrest = params.ncside % params.xSize;
	params.sizeVertical = params.ncside / params.ySize;
	long sizeVerticalrest = params.ncside % params.ySize;

	/*if(rank == 0) {
		printf("nr proc: %d, x: %ld, y: %ld, x*y = %ld\nsizeVertical: %ld sizeVerticalrest: %ld, sizeHorizontal: %ld sizeHorizontalrest: %ld\n\n", numberOfProcess, params.xSize, params.ySize,  params.xSize * params.ySize, params.sizeVertical, sizeVerticalrest, params.sizeHorizontal, sizeHorizontalrest);fflush(stdout);
	}*/

	params.sizeHorizontalBase = params.sizeHorizontal + (sizeHorizontalrest > 0);
	params.sizeVerticalBase = params.sizeVertical + (sizeVerticalrest > 0);

	if((rank % params.xSize)  < sizeHorizontalrest) {
		params.sizeHorizontal++;
		params.xLowerBound = (rank % params.xSize) * params.sizeHorizontal;
	}
	else {
		params.xLowerBound = (rank % params.xSize) * params.sizeHorizontal + sizeHorizontalrest;
	}
	

	if((rank / params.xSize) < sizeVerticalrest) {
		params.sizeVertical++;
		params.yLowerBound = (rank / params.xSize ) * params.sizeVertical;
	}
	else {
		params.yLowerBound = (rank / params.xSize ) * params.sizeVertical + sizeVerticalrest;
	}

	params.xUpperBound = params.xLowerBound + params.sizeHorizontal - 1;
	
	params.yUpperBound = params.yLowerBound + params.sizeVertical - 1;

	//printf("rank %d\nxmin: %ld, xmax: %ld\nymin: %ld, ymax: %ld\nsizeVertical: %ld, sizeHorizontal: %ld\n\n", rank, params.xLowerBound, params.xUpperBound, params.yLowerBound, params.yUpperBound, params.sizeVertical, params.sizeHorizontal);fflush(stdout);

	return numberOfProcess;
}

grid_t initTotalGrid(grid_t grid, long ncside) {
	grid.centerOfMassX = (double*) malloc(ncside*ncside*sizeof(double));
	if(grid.centerOfMassX == NULL) {
		printf("ERROR malloc centerOfMassX\n");fflush(stdout);
		exit(0);
	}
	grid.centerOfMassY = (double*) malloc(ncside*ncside*sizeof(double));
	if(grid.centerOfMassY == NULL) {
		printf("ERROR malloc centerOfMassY\n");fflush(stdout);
		exit(0);
	}

	grid.m = (double *) malloc(ncside*ncside*sizeof(double));
	if(grid.m == NULL) {
		printf("ERROR malloc m\n");fflush(stdout);
		exit(0);
	}
	
	return grid;
}

grid_tt ** initGridSendReceive(int rank) {
	grid_tt **grid;
	grid = (grid_tt **) malloc(16*sizeof(grid_tt *));
	if(grid == NULL) {
		printf("ERROR malloc GridSendReceive \n");fflush(stdout);
		exit(0);
	}
	for (int i = 0; i < 16; i += 4) {
		grid[i] = (grid_tt *) malloc(params.sizeVertical*sizeof(grid_tt));
		if(grid[i] == NULL) {
			printf("ERROR malloc %d-1 \n", i);fflush(stdout);
			exit(0);
		}
		grid[i+1] = (grid_tt *) malloc(sizeof(grid_tt));
		if(grid[i+1] == NULL) {
			printf("ERROR malloc %d-2 \n", i);fflush(stdout);
			exit(0);
		}
		grid[i+2] = (grid_tt *) malloc(params.sizeHorizontal*sizeof(grid_tt));
		if(grid[i+2] == NULL) {
			printf("ERROR malloc %d-3 \n", i);fflush(stdout);
			exit(0);
		}
		grid[i+3] = (grid_tt *) malloc(sizeof(grid_tt));
		if(grid[i+3] == NULL) {
			printf("ERROR malloc %d-4 \n", i);fflush(stdout);
			exit(0);
		}
	}

	return grid;
}

int * findNeighborsRank(int* idToSend, int rank, int numberOfProcess) {
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

	return(idToSend);
}

void freeEverything(particle_t *par, grid_t particleGrid, long long nside){
	free(par);

	free(particleGrid.m);
	free(particleGrid.centerOfMassX);
	free(particleGrid.centerOfMassY);
	return;
}

particle_t_reduced * initParReceived(long long n_part, long *size) {
	particle_t_reduced *par;

	if(n_part <= 100) *size = n_part;
	else if(n_part <= 500) *size = n_part/2;
	else *size = n_part/ 10;
	par = (particle_t_reduced *)malloc((*size)*sizeof(particle_t_reduced));
	if(par == NULL) {
		printf("ERROR malloc Par Received\n");fflush(stdout);
		exit(0);
	}

	return(par);
}	

particle_t_reduced ** initParSend(long long n_part, long *size) {
	particle_t_reduced **par;

	if(n_part <= 100) *size = n_part;
	else if(n_part <= 500) *size = n_part/2;
	else *size = n_part/ 10;

	par = (particle_t_reduced **) malloc(8*sizeof(particle_t_reduced *));
	if(par == NULL) {
		printf("ERROR malloc GridSendReceive \n");fflush(stdout);
		exit(0);
	}
	for (int i = 0; i < 8; i++) {
		par[i] = (particle_t_reduced *) malloc((*size)*sizeof(particle_t_reduced));
		if(par[i] == NULL) {
			printf("ERROR malloc ParSend %d\n", i);fflush(stdout);
			exit(0);
		}
	}

	return(par);
}	

