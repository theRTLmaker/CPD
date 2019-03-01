#ifndef INIT_PROGRAM_H
#define INIT_PROGRAM_H
#include "physics.h"


particle_t * handler_input(int argc ,char *argv[]);

void init_particles(long seed, long ncside, long long n_part, particle_t *par);

particle_t * CreateParticleArray(long long int n_part);

particle_t *** CreateParticleGrid(long long int n_part);

grid_t initGrid(grid_t grid, long ncside);

#endif