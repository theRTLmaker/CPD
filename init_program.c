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

particle_t * CreateParticleArray(long long n_part) {
	particle_t *par = (particle_t*) malloc(n_part*sizeof(particle_t));
	if(par ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	return par;
}

void init_particles(particle_t *par) {
    long long i;

    srandom(params.seed);

    for(i = 0; i < params.n_part; i++) {        
		par[i].positionX = RND0_1;
        par[i].positionY = RND0_1;
        
        par[i].vx = RND0_1 / params.ncside / 10.0;
        par[i].vy = RND0_1 / params.ncside / 10.0;

        par[i].m = RND0_1 * params.ncside / (G * 1e6 * params.n_part);

        par[i].gridCoordinateX = par[i].positionX * params.ncside / 1;
        par[i].gridCoordinateY = par[i].positionY * params.ncside / 1;

        par[i].number = -1;
    }
}

// Talvez para tentar paralelizar
void computeGridPosition(particle_t *par) {
	long long i;

	for(i = 0; i < params.n_part; i++) {        
		par[i].gridCoordinateX = par[i].positionX * params.ncside / 1;
        par[i].gridCoordinateY = par[i].positionY * params.ncside / 1;

        par[i].number = i;
    }
}

int findGridDivision(int numberOfProcess, int rank) {
	int dims[2] = {0,0};
	MPI_Dims_create(numberOfProcess, 2, dims);
	if(dims[1] == 1) {
		do {
			numberOfProcess--;
			MPI_Dims_create(numberOfProcess, 2, dims);
		}while(dims[1] == 1);
	}
	params.xSize = dims[1];
	params.ySize = dims[0];

	if(rank == 0)
		printf("nr proc: %d, x: %ld, y: %ld, x*y = %ld\n", numberOfProcess, params.xSize, params.ySize,  params.xSize * params.ySize);


	// Divisao em grelha pelos processos
	if(rank != numberOfProcess - 1) {
		// Limite inferior do X
		params.xLowerBound = (rank % params.xSize) * params.ncside/params.xSize;

		// Limite inferior do Y
		params.yLowerBound = (rank / params.ySize) * params.ncside/params.ySize;

		// Quando o processo esta no limite X da grelha
		if(rank % params.xSize == params.xSize - 1) {
			// limite superior do X 
			params.xUpperBound = params.ncside - 1;

			// Quando o processo esta no limite Y da grelha
			if(rank / params.ySize == params.ySize - 1) 
				params.yUpperBound = params.ncside - 1;
			else
				params.yUpperBound = ((rank + 1) / params.ySize) * params.ncside/params.ySize - 1;
		}
		else{
			params.xUpperBound = ((rank + 1) % params.xSize) * params.ncside/params.xSize - 1;

			if(rank / params.ySize == params.ySize - 1)
				params.yUpperBound = params.ncside - 1;
			else
				params.yUpperBound = ((rank + 1) / params.ySize) * params.ncside/params.ySize - 1;
		}
	}

	params.sizeSmallD = params.xUpperBound - params.xLowerBound;
	params.sizeBigD = params.yUpperBound - params.yLowerBound;

	printf("rank %d\nxmin: %ld, xmax: %ld\nymin: %ld, ymax: %ld\n", rank, params.xLowerBound, params.xUpperBound, params.yLowerBound, params.yUpperBound);
	return numberOfProcess;
}

grid_t initTotalGrid(grid_t grid, long ncside) {
	grid.centerOfMassX = (double*) malloc(ncside*ncside*sizeof(double));
	if(grid.centerOfMassX == NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	grid.centerOfMassY = (double*) malloc(ncside*ncside*sizeof(double));
	if(grid.centerOfMassY == NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	grid.m = (double *) malloc(ncside*ncside*sizeof(double));
	if(grid.m == NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	
	return grid;
}

grid_tt ** initGridSendReceive(int rank) {
	grid_tt **grid;
	grid = (grid_tt **) malloc(8*sizeof(grid_tt *));
	if(grid == NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	for (int i = 0; i < 4; ++i) {
		grid[i] = (grid_tt *) malloc((params.yUpperBound-params.yLowerBound + 1)*sizeof(grid_tt));
		if(grid[i] == NULL) {
			printf("ERROR malloc\n");
			exit(0);
		}
		grid[i+1] = (grid_tt *) malloc(sizeof(grid_tt));
		if(grid[i+1] == NULL) {
			printf("ERROR malloc\n");
			exit(0);
		}
		grid[i+2] = (grid_tt *) malloc((params.xUpperBound-params.xLowerBound + 1)*sizeof(grid_tt));
		if(grid[i+2] == NULL) {
			printf("ERROR malloc\n");
			exit(0);
		}
		grid[i+3] = (grid_tt *) malloc(sizeof(grid_tt));
		if(grid[i+3] == NULL) {
			printf("ERROR malloc\n");
			exit(0);
		}
	}

	return grid;
}

/*grid_t initPartialGrid(int numberOfProcess, int processID, grid_t grid, parameters *params) {
	// Divisão em linhas pelos processos
	if(params.xSize == 1) {
		if(processID == numberOfProcess - 1)
			params.gridSize = params.ncside*(params.ncside/numberOfProcess + params.ncside%numberOfProcess);
		else
			params.gridSize = params.ncside*params.ncside/numberOfProcess;
	}
	// Divisao em grelha pelos processos
	else {
		if(processID % params.xSize == params.xSize - 1) {
			if(processID % params.ySize == params.ySize - 1) 
				params.gridSize = ((params.ncside/params.xSize) + (params.ncside%params.xSize)) * ((params.ncside/params.ySize) + (params.ncside%params.ySize));
			else
				params.gridSize = ((params.ncside/params.xSize) + (params.ncside%params.xSize)) * (params.ncside/params.ySize);
		}
		else{
			if(processID % params.ySize == params.ySize - 1) 
				params.gridSize = (params.ncside/params.xSize) * ((params.ncside/params.ySize) + (params.ncside%params.ySize));
			else
				params.gridSize = (params.ncside/params.xSize)*(params.ncside/params.ySize);
		}
	}

	grid.centerOfMassX = (double*) malloc(params.gridSize*sizeof(double));
	if(grid.centerOfMassX ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	grid.centerOfMassY = (double*) malloc(params.gridSize*sizeof(double));
	if(grid.centerOfMassY ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	grid.m = (double *) malloc(params.gridSize*sizeof(double));
	if(grid.m ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	
	return grid;
}
*/

void freeEverything(particle_t *par, grid_t particleGrid, long long nside){
	free(par);

	free(particleGrid.m);
	free(particleGrid.centerOfMassX);
	free(particleGrid.centerOfMassY);
	return;
}
