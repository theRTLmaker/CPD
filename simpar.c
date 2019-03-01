#include <stdio.h>
#include <stdlib.h>
#include "physics.h"
#include "init_program.h"

/* Fazer um handler para tratar dos argumentos de entrada e passar para o init_particles*/

int main(int argc, char *argv[])
{
	grid_t grid;
	particle_t *par;
	particle_t center;
	int nr_part; 
	int i = 0;
	int j = 0;

	par = handler_input(argc, argv);
	nr_part = atoll(argv[3]);
	long long int timeStep = atoll(argv[4]);
	



	//Tá com um erro esquisito
	for (j = 0; j < timeStep; j++)
	{
		center = calculateCenterOfMass(par, nr_part); 

		for(i = 0; i<nr_part; i++){
			par[i].appliedForce = calculateGravForce(par[i], center);
			par[i].position = calculateNextPosition(par[i]);
			par[i].velocity = calculateNextVelocity(par[i]);
		}

		printAllParticles(par, nr_part);
		
		printCenter(center);
	}

	freeEverything(par, grid.cells, nr_part);
	return 0;
}

