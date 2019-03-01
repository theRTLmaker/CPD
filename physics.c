#include "physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double vectorNorm(vector2 a){
	return sqrt(a.x * a.x + a.y * a.y);
}

vector2 addVectors(vector2 a, vector2 b){
	vector2 r;
	r.x = a.x + b.x;
	r.y = a.y + b.y;
	return r;
}

vector2 subVectors(vector2 a, vector2 b){
	vector2 r;
	r.x = a.x - b.x;
	r.y = a.y - b.y;
	return r;
}

vector2 multiplyVectorByConst(double c, vector2 v){
	vector2 r;
	r.x = c*v.x;
	r.y = c*v.y;
	return r;
}

particle_t calculateCenterOfMass(particle_t *par, long long n_part){
	long long i;
	particle_t center;
	center.position.x = 0;
	center.position.y = 0;
	center.m = 0;
	for(i = 0; i < n_part; i++)
    {
        center.position = addVectors(center.position, multiplyVectorByConst(par[i].m, par[i].position));
        center.m += par[i].m;
    }
    center.position = multiplyVectorByConst(1/center.m, center.position);
    return center;
}

vector2 calculateGravForce(particle_t p1, particle_t massCenter){
	vector2 gravForce;
	double gravForceMag;

	vector2 forceDirection = subVectors(massCenter.position, p1.position);

	double distance = vectorNorm(forceDirection);
	if(distance < EPSLON)
		gravForceMag = 0;
	else
		gravForceMag = (p1.m * massCenter.m * G) /(distance*distance);

	return multiplyVectorByConst(gravForceMag / distance, forceDirection);
}

vector2 calculateNextPosition(particle_t particle){ // x = x0 + v0t + 0.5 a t^2 (t = 1) a = F/m
	vector2 a = multiplyVectorByConst(1/particle.m, particle.appliedForce); //a = F/m
	vector2 newPos = addVectors(particle.position, addVectors(particle.velocity, multiplyVectorByConst(0.5, a))); //x = x0 + v0t + 0.5 a t^2 (t = 1)

	if(newPos.x >= 1) newPos.x = newPos.x - floor(newPos.x);
	else if(newPos.x < 0) newPos.x = 1 + (newPos.x - ceil(newPos.x)); 

	if(newPos.y >= 1) newPos.y = newPos.y - floor(newPos.x);
	else if(newPos.y < 0) newPos.y = 1 + (newPos.y - ceil(newPos.y));

	return newPos;
} 

vector2 calculateNextVelocity(particle_t particle){ // v = v0 + at (t = 1)
	vector2 a = multiplyVectorByConst(1/particle.m, particle.appliedForce); //a = F/m
	vector2 newVel = addVectors(particle.velocity, a);
	return newVel;
}

void printVectorPosition(vector2 p){
	printf("	Position - ( %lf , %lf )\n", p.x, p.y);
}

void printVectorVelocity(vector2 v) {
	printf("	Velocity - ( %lf , %lf )\n",v.x, v.y);
}

void printVectorGrid(vector2grid g) {
	printf("	Grid - ( %d , %d )\n",g.x, g.y);
}



void printParticle(particle_t p) {
	printVectorPosition(p.position);
	printVectorVelocity(p.velocity);
	printVectorGrid(p.gridCoordinate);
	printf("	Force - ( %lf , %lf )\n", p.appliedForce.x, p.appliedForce.y);
	printf("	mass: %1f\n", p.m);
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
	printf("	mass: %1f\n", p.m);
}