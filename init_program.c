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

particle_t * CreateParticleArray(long long int n_part) {
	particle_t *par = (particle_t*) malloc(n_part*sizeof(particle_t));
	if(par ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	return par;
}

grid_t initGrid(grid_t grid, long ncside) {
	grid.centerOfMass = (vector2**) malloc(ncside*sizeof(vector2 *));
	if(grid.centerOfMass ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	grid.m = (long double **) malloc(ncside*sizeof(long double *));
	if(grid.m ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}
	for(int i = 0; i < ncside; i++) {
		grid.centerOfMass[i] = (vector2 *) malloc(ncside*sizeof(vector2));
		if(grid.centerOfMass[i] == NULL) {
			printf("ERROR malloc\n");
			exit(0);
		}
		grid.m[i] = (long double *) malloc(ncside*sizeof(long double));
		if(grid.m[i] == NULL) {
			printf("ERROR malloc\n");
			exit(0);
		}
	}
	
	return grid;
}

vector2grid findPosition(particle_t par, long ncside) {
	vector2grid vector;
	vector.x = par.position.x * ncside / 1;
	vector.y = par.position.y * ncside / 1;

	return vector;
}


void init_particles(long seed, long ncside, long long int n_part, particle_t *par) {
    long long i;

    srandom(seed);

    for(i = 0; i < n_part; i++) {        
		par[i].position.x = RND0_1;
        par[i].position.y = RND0_1;
        
        par[i].velocity.x = RND0_1 / ncside / 10.0;
        par[i].velocity.y = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);

        par[i].gridCoordinate = findPosition(par[i], ncside);

        par[i].pastPositions = (vector2 *)malloc(1000* sizeof(vector2));
    }
}


void freeEverything(particle_t *par, grid_t particleGrid, long long int nside){
	free(par);

	for(int i = 0; i < nside; i++) {
			free(particleGrid.m[i]);
			free(particleGrid.centerOfMass[i]);
	}

	free(particleGrid.m);
	free(particleGrid.centerOfMass);

	return;
}
