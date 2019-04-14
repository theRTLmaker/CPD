#include <stdio.h>
#include "physics.h"
#include "debug.h"

/*void printVectorPosition(vector2 p){
	printf("	Position - ( %f , %f )\n", p.x, p.y);
}

void printVectorVelocity(vector2 v) {
	printf("	Velocity - ( %f , %f )\n",v.x, v.y);
}

void printVectorGrid(vector2grid g) {
	printf("	Grid - ( %d , %d )\n",g.x, g.y);
}*/

void printParticle(particle_t p) {
	//printVectorGrid(p.gridCoordinate);
	printf("	Position - ( %f , %f )\n", p.positionX, p.positionY);
	printf("	Velocity - ( %f , %f )\n", p.vx, p.vy);
	printf("	Force - ( %f , %f )\n", p.appliedForceX, p.appliedForceY);
	printf("	mass: %f\n", p.m);
}

void printAllParticles(particle_t *p, long long int nr_part) {
	long long int i = 0;
	while(i < nr_part) 
	{	
		printf("Particle - %lld\n", i);
		printParticle(p[i]);
		i++;
	}
}

void printCenter(particle_t p) {
	printf("Center of Mass\n");
	printf("	Position - ( %f , %f )\n", p.positionX, p.positionY);
	printf("	Velocity - ( %f , %f )\n", p.vx, p.vy);
	printf("	mass: %f\n", p.m);
}
