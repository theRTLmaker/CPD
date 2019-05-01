#include "physics.h"
#include <math.h>


void calculateGravForce(particle_t *p1, double massCenterX, double massCenterY, float m, int sideUPDOWN, int sideLEFTRIGHT){
	double gravForceMag;

	double forceDirectionX;
	double forceDirectionY;
	double invDist;
	switch(sideUPDOWN) {
		case(UP):
			massCenterY = massCenterY + 1;
			break;
		case(DOWN):
			massCenterY = massCenterY - 1;
			break;
		default:
			break;
	}

	switch(sideLEFTRIGHT) {
		case(LEFT):
			massCenterX = massCenterX - 1;
			break;
		case(RIGHT):
			massCenterX = massCenterX + 1;
			break;
		default:
			break;
	}
	//massCenterX = massCenterX + sideUPDOWN;
	//massCenterY = massCenterY + sideLEFTRIGHT;

	forceDirectionX = massCenterX - p1 -> positionX;
	forceDirectionY = massCenterY - p1 -> positionY;

	double distance = vectorNorm(forceDirectionX, forceDirectionY);
	if(distance < EPSLON)
		return;
	else{
		invDist = 1.0/distance;
		gravForceMag = (p1 -> m * m * G) * invDist * invDist;
		p1 -> appliedForceX += gravForceMag * forceDirectionX * invDist;
		p1 -> appliedForceY += gravForceMag * forceDirectionY * invDist;
	}

	return;
}


