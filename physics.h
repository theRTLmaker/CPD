#ifndef PHYSICS_H
#define PHYSICS_H

#include <stdio.h>
#include <stdlib.h>

#define G 6.67408e-11
#define EPSLON 0.0005

#define MIDDLE -1
#define DOWN 0
#define UP 1
#define LEFT 2
#define RIGHT 3


typedef struct _vector2{
	float x;
	float y;
} vector2;

typedef struct _vector2grid{
	int x;
	int y;
} vector2grid;

typedef struct _particle_t {
	float m;
	float positionX;
	float positionY;
	vector2 velocity;
	int gridCoordinateX;
	int gridCoordinateY;
	vector2 appliedForce;
	vector2 *pastPositions;
} particle_t;

typedef struct _grid_t{
	float **m;
	vector2 **centerOfMass;
	unsigned char **mask;
} grid_t;

#define vectorNorm(x, y) (sqrt(x * x + y * y))
#define constrain1(size, n) ((n) == (-1) ? ((n) + (size)) : (n))
#define constrain2(size, n) ((n) == (size) ? ((n) - (size)) : (n))
#define SUM_A( x, y )  ((x) == 0 || (y) == 0 ? 0 : ( ( ( (x) * (x) ) / ( ( x ) + ( y ) ) ) * ( y ) ))
#define compareVectorsGrid(ax, ay, bx, by) ((ax) == (bx) && (ay) == (by) ? 1 : 0)

//int constrain(int size, int n);
float addVectors(float a, float b);
vector2 subVectors(vector2 a, vector2 b);
vector2 multiplyVectorByConst(float c, vector2 v);

vector2 calculateGravForce(particle_t p1, vector2 massCenter, float m, int sideUPDOWN, int sideLEFTRIGHT);

vector2 calculateNextPosition(particle_t particle); // x = x0 + v0t + 0.5 a t^2 (t = 1)

vector2 calculateNextVelocity(particle_t particle); // v = v0 + at (t = 1)


#endif