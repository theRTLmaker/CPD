#include "init_program.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND0_1 ((double) random() / ((long long)1<<31))



parameters handler_input(int argc ,char *argv[]) 
{
	parameters params;
	long seed = 0;

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

	if(params.seed < 0 || params.ncside < 0 || params.n_part < 0){
		printf("Wrong parameteres\n");
		exit(0);
	}

	return params;
}

particle_t * CreateParticleArray(long long n_part) {
	particle_t *par = (particle_t*) malloc(n_part*sizeof(particle_t));
	if(par ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	return par;
}

void init_particles(parameters params, particle_t *par) {
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

        par[i].active = 1;
    }
}

// Talvez para tentar paralelizar
void computeGridPosition(parameters params, particle_t *par) {
	long long i;

	for(i = 0; i < params.n_part; i++) {        
		par[i].gridCoordinateX = par[i].positionX * params.ncside / 1;
        par[i].gridCoordinateY = par[i].positionY * params.ncside / 1;

        par[i].active = 1;
    }
}

void findGridDivision(parameters *params, int numberOfProcess) {
	long aux;
	if(sqrt(numberOfProcess) % 1.0 == 0) {
		params->xSize = sqrt(numberOfProcess);
		params->ySize = sqrt(numberOfProcess);
	}
	else {
		params->xSize = sqrt(numberOfProcess)/1;
		params->ySize = params->xSize;
		aux = params->xSize * params->ySize;
		do {
			if(aux < numberOfProcess) {
				params->ySize += 1;
			}
			else if(aux > numberOfProcess){
				params->xSize -= 1;
			}
			aux = params->xSize * params->ySize;

		}while (aux != numberOfProcess)
	}
	printf("nr proc: %d, x %ld, y %ld\n, x*y %ld", numberOfProcess, params->xSize, params->ySize, aux);
}

grid_t initTotalGrid(grid_t grid, long ncside) {
	grid.centerOfMassX = (double*) malloc(ncside*ncside*sizeof(double ));
	if(grid.centerOfMassX ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	grid.centerOfMassY = (double*) malloc(ncside*ncside*sizeof(double ));
	if(grid.centerOfMassY ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	grid.m = (double *) malloc(ncside*ncside*sizeof(double ));
	if(grid.m ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	
	return grid;
}

grid_t initPartialGrid(int numberOfProcess, int processID, grid_t grid, long ncside) {
	grid.centerOfMassX = (double*) malloc(ncside*ncside*sizeof(double ));
	if(grid.centerOfMassX ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	grid.centerOfMassY = (double*) malloc(ncside*ncside*sizeof(double ));
	if(grid.centerOfMassY ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	grid.m = (double *) malloc(ncside*ncside*sizeof(double ));
	if(grid.m ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	
	return grid;
}


void freeEverything(particle_t *par, grid_t particleGrid, long long nside){
	free(par);

	free(particleGrid.m);
	free(particleGrid.centerOfMassX);
	free(particleGrid.centerOfMassY);
	return;
}
