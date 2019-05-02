#ifndef INIT_PROGRAM_H
#define INIT_PROGRAM_H
#include "physics.h"

#define LEFTPROCESS 0
#define	UPLEFTPROCESS 1
#define	UPPROCESS 2
#define	UPRIGHTPROCESS 3
#define	RIGHTPROCESS 4
#define	DOWNRIGHTPROCESS 5
#define	DOWNPROCESS 6
#define	DOWNLEFTPROCESS 7

typedef struct _parameters{
	long seed;
	long ncside;
	long long gridSize;
	long long n_part;
	long timeStep;
	long xSize;
	long ySize;
	long xLowerBound;
	long xUpperBound;
	long yLowerBound;
	long yUpperBound;
	long sizeBigD;
	long sizeSmallD;
} parameters;

extern parameters params;

void handler_input(int argc ,char *argv[]);

vector2grid findPosition(particle_t par, long ncside);

void init_particles(particle_t *par);

particle_t * CreateParticleArray(long long n_part);

// Talvez para tentar paralelizar
void computeGridPosition(particle_t *par);

void findGridDivision(int numberOfProcess, int rank);

grid_tt ** initGridSendReceive(int rank);

grid_t initTotalGrid(grid_t grid, long ncside);

grid_t initPartialGrid(int numberOfProcess, int processID, grid_t grid);

void freeEverything(particle_t *par, grid_t particleGrid, long long int nside);

#endif