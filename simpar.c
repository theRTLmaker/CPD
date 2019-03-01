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

	par = handler_input(argc, argv);
	nr_part = atoll(argv[3]);

	printAllParticles(par, nr_part);

	center = calculateCenterOfMass(par, nr_part); 

	printCenter(center);
	return 0;
}

