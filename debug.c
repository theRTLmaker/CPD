#include <stdio.h>
#include "physics.h"
#include "debug.h"

void printVectorPosition(vector2 p){
	printf("	Position - ( %f , %f )\n", p.x, p.y);
}

void printVectorVelocity(vector2 v) {
	printf("	Velocity - ( %f , %f )\n",v.x, v.y);
}

void printVectorGrid(vector2grid g) {
	printf("	Grid - ( %d , %d )\n",g.x, g.y);
}



void printParticle(particle_t p) {
	printVectorPosition(p.position);
	printVectorVelocity(p.velocity);
	//printVectorGrid(p.gridCoordinate);
	printf("	Force - ( %f , %f )\n", p.appliedForce.x, p.appliedForce.y);
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
	printVectorPosition(p.position);
	printVectorVelocity(p.velocity);
	printf("	mass: %f\n", p.m);
}
