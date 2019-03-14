#include "physics.h"
#include <math.h>

inline vector2 addVectors(vector2 a, vector2 b) {
	vector2 r;
	r.x = a.x + b.x;
	r.y = a.y + b.y;
	return r;
}

inline vector2 subVectors(vector2 a, vector2 b){
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

vector2 calculateGravForce(particle_t p1, vector2 massCenter, long double m, int sideUPDOWN, int sideLEFTRIGHT){
	double gravForceMag;

	vector2 forceDirection;

	if(m == 0) {
		forceDirection.x = 0;
		forceDirection.y = 0;
		return forceDirection;
	}
	vector2 aux = massCenter;

	switch(sideUPDOWN) {
		case(UP):
			aux.y += 1;
			break;
		case(DOWN):
			aux.y -= 1;
			break;
		default:
			break;
	}

	switch(sideLEFTRIGHT) {
		case(LEFT):
			aux.x -= 1;
			break;
		case(RIGHT):
			aux.x += 1;
			break;
		default:
			break;
	}

	forceDirection = subVectors(aux, p1.position);

	double distance = vectorNorm(forceDirection);
	if(distance < EPSLON){
		gravForceMag = 0;
		distance = 1;
	}
	else
		gravForceMag = (p1.m * m * G) /(distance*distance);

	return multiplyVectorByConst(gravForceMag / distance, forceDirection);
}

vector2 calculateNextPosition(particle_t particle){ // x = x0 + v0t + 0.5 a t^2 (t = 1) a = F/m
	vector2 a = multiplyVectorByConst(1/particle.m, particle.appliedForce); //a = F/m
	vector2 newPos = addVectors(particle.position, addVectors(particle.velocity, multiplyVectorByConst(0.5, a))); //x = x0 + v0t + 0.5 a t^2 (t = 1)

	if(newPos.x >= 1) newPos.x = newPos.x - floor(newPos.x);
	else if(newPos.x < 0) newPos.x = 1 + (newPos.x - ceil(newPos.x)); 

	if(newPos.y >= 1) newPos.y = newPos.y - floor(newPos.y);
	else if(newPos.y < 0) newPos.y = 1 + (newPos.y - ceil(newPos.y));

	return newPos;
} 

vector2 calculateNextVelocity(particle_t particle){ // v = v0 + at (t = 1)
	vector2 a = multiplyVectorByConst(1/particle.m, particle.appliedForce); //a = F/m
	vector2 newVel = addVectors(particle.velocity, a);
	return newVel;
}

