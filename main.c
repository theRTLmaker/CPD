#include <stdio.h>
#include <stdlib.h>
#include "physics.h"

int main(int argc, char const *argv[])
{
	int n = 2;
	particle_t par[2];


	//init_particles(30000, 10, )

	particle_t a, b, c;
	a.position.x = 1;
	a.position.y = 2;
	b.position.x = 5;
	b.position.y = 6;
	a.m = 2;
	b.m = 5;
	par[0] = a;
	par[1] = b;

	c = calculateCenterOfMass(par, n);
	printParticle(a);
	printParticle(b);
	printParticle(c);

	return 0;
}