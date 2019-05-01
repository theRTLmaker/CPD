#ifndef INIT_PROGRAM_H
#define INIT_PROGRAM_H
#include "physics.h"



typedef struct _parameters{
	long seed;
	long ncside;
	long long gridSize;
	long long n_part;
	long long partialNrPart;
	long timeStep;
	long xSize;
	long ySize;
} parameters;

parameters handler_input(int argc ,char *argv[]);

vector2grid findPosition(particle_t par, long ncside);

void init_particles(parameters params, particle_t *par);

particle_t * CreateParticleArray(long long n_part);

// Talvez para tentar paralelizar
void computeGridPosition(parameters params, particle_t *par);

void findGridDivision(parameters *params, int numberOfProcess);

grid_t initTotalGrid(grid_t grid, long ncside);

grid_t initPartialGrid(int numberOfProcess, int processID, grid_t grid, parameters *params);

void freeEverything(particle_t *par, grid_t particleGrid, long long int nside);

#endif