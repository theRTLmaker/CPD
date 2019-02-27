#include "physics.h"
#include <stdio.h>
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
	gravForceMag = (p1.m * massCenter.m * G) /(distance*distance);
	return multiplyVectorByConst(gravForceMag / distance, forceDirection);
}


void printVector(vector2 v);

void printParticle(particle_t p){
	printf("p: ");
	printVector(p.position);
	printf("v: ");
	printVector(p.velocity);	
	printf("m: %f\n", p.m);
}



void printVector(vector2 v){
	printf("(%f, %f)\n", v.x, v.y);
}