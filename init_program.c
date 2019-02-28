#include "init_program.h"
#include <stdio.h>
#include <stdlib.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define EPSLON 0.01

particle_t *CreateParticleArray(long long int n_part)
{
	particle_t *par = (particle_t*) malloc(n_part*sizeof(particle_t));
	if(par ==NULL) {
		printf("ERROR malloc\n");
		exit(0);
	}

	return par;
}


void ParameteresControl (int argc, long seed, long ncside, long long int n_part)
{
	if(argc != 5){
		printf("Wrong parameteres\n");
		exit(0);
	}
	else if(seed < 0 || ncside < 0 || n_part < 0){
		printf("Wrong parameteres\n");
		exit(0);
	}

}

particle_t *handler_input(int argc ,char *argv[])
{
	particle_t *par;
	long seed = 0;
	long ncside = 0;
	long long int n_part = 0;

	seed = atol(argv[1]);
	//printf("seed = %ld\n", seed );
	ncside = atol(argv[2]);
	//printf("ncside = %ld\n", ncside );
	n_part = atoll(argv[3]);
	//printf("n_part = %lld\n", n_part);
	//time step ?????

	ParameteresControl(argc, seed, ncside, n_part);

	par = CreateParticleArray(n_part);

	init_particles(seed, ncside, n_part, par);

	return par;

}

void init_particles(long seed, long ncside, long long int n_part, particle_t *par)
{
    long long i;

    srandom(seed);

    for(i = 0; i < n_part; i++)
    {
        par[i].position.x = RND0_1;
        par[i].position.y = RND0_1;
        par[i].velocity.x = RND0_1 / ncside / 10.0;
        par[i].velocity.y = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}
