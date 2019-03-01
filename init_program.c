#include "init_program.h"
#include <stdio.h>
#include <stdlib.h>

#define RND0_1 ((double) random() / ((long long)1<<31))



particle_t *handler_input(int argc ,char *argv[]) {
	particle_t *par;
	long seed = 0;
	long ncside = 0;
	long long int n_part = 0;
	long long int timeStep = 0;

	if(argc != 5){
		printf("Wrong number of input parameteres\n");
		exit(0);
	}

	seed = atol(argv[1]);
	//printf("seed = %ld\n", seed );
	ncside = atol(argv[2]);
	//printf("ncside = %ld\n", ncside );
	n_part = atoll(argv[3]);
	//printf("n_part = %lld\n", n_part);
	//time step ?????
	timeStep = atoll(argv[4]);

	if(seed < 0 || ncside < 0 || n_part < 0){
		printf("Wrong parameteres\n");
		exit(0);
	}

	par = CreateParticleArray(n_part);

	init_particles(seed, ncside, n_part, par);

	return par;
}

particle_t * CreateParticleArray(long long int n_part) {
	particle_t *par = (particle_t*) malloc(n_part*sizeof(particle_t));
	if(par ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	//par->nextParticle = NULL;

	return par;
}


particle_t *** CreateParticleGrid(long long int n_part) {
	particle_t ***cells = (particle_t***) malloc(n_part*sizeof(particle_t**));
	if(cells ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	for(int i = 0; i < n_part; i++) {
		cells[i] = (particle_t**) malloc(n_part*sizeof(particle_t*));
		if(cells[i] == NULL) {
			printf("ERROR malloc\n");
			exit(0);
		}
	}

	return cells;
}

grid_t initGrid(grid_t grid, long ncside) {
	grid.dimension = ncside;
	grid.cells = CreateParticleGrid(ncside);

	return grid;
}

particle_t findPosition(particle_t par, long ncside) {
	par.gridCoordinate.x = par.position.x * ncside / 1;
	par.gridCoordinate.y = par.position.y * ncside / 1;

	return par;
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

        par[i] = findPosition(par[i], ncside);
    }
}


void freeEverything(particle_t *par, particle_t ***particleGrid, long long int n_part){
	free(par);

	/*for(int i = 0; i < n_part; i++) {
		free(particleGrid[i]);
	}*/

	//free(particleGrid);

	return;
}
