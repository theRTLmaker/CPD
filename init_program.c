#include "init_program.h"
#include <stdio.h>
#include <stdlib.h>

#define RND0_1 ((double) random() / ((long long)1<<31))



particle_t *handler_input(int argc ,char *argv[], parameters *param) {
	particle_t *par;
	long seed = 0;

	if(argc != 5){
		printf("Wrong number of input parameteres\n");
		exit(0);
	}

	seed = atol(argv[1]);
	//printf("seed = %ld\n", seed );
	param -> ncside = atol(argv[2]);
	//printf("ncside = %ld\n", ncside );
	param -> n_part = atoll(argv[3]);
	//printf("n_part = %lld\n", n_part);
	//time step ?????
	param -> timeStep = atoll(argv[4]);

	if(seed < 0 || param -> ncside < 0 || param -> n_part < 0){
		printf("Wrong parameteres\n");
		exit(0);
	}

	par = CreateParticleArray(param -> n_part);

	init_particles(seed, param -> ncside, param -> n_part, par);

	return par;
}

particle_t * CreateParticleArray(long long n_part) {
	particle_t *par = (particle_t*) malloc(n_part*sizeof(particle_t));
	if(par ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	return par;
}

grid_t initGrid(grid_t grid, long ncside) {
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

/*vector2grid findPosition(particle_t par, long ncside) {
	vector2grid vector;
	vector.x = par.position.x * ncside / 1;
	vector.y = par.position.y * ncside / 1;

	return vector;
}*/


void init_particles(long seed, long ncside, long long int n_part, particle_t *par) {
    long long i;

    srandom(seed);

    for(i = 0; i < n_part; i++) {        
		par[i].positionX = RND0_1;
        par[i].positionY = RND0_1;
        
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);

        par[i].gridCoordinateX = par[i].positionX * ncside / 1;
        par[i].gridCoordinateY = par[i].positionY * ncside / 1;

    }
}


void freeEverything(particle_t *par, grid_t particleGrid, long long nside){
	free(par);

	free(particleGrid.m);
	free(particleGrid.centerOfMassX);
	free(particleGrid.centerOfMassY);
	return;
}
