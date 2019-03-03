#ifndef INIT_PROGRAM_H
#define INIT_PROGRAM_H
#include "physics.h"



typedef struct _parameters{
	long ncside;
	long long n_part;
	long long timeStep;
} parameters;

particle_t * handler_input(int argc ,char *argv[], parameters *param);

void init_particles(long seed, long ncside, long long n_part, particle_t *par);

particle_t * CreateParticleArray(long long int n_part);

gridCell ** CreateParticleGrid(long long int n_part);

grid_t initGrid(grid_t grid, long ncside);

void freeEverything(particle_t *par, gridCell **particleGrid, long long int nside);

#endif