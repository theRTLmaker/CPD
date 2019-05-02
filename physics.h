#ifndef PHYSICS_H
#define PHYSICS_H

#include <stdio.h>
#include <stdlib.h>

#define G 6.67408e-11
#define EPSLON 0.0005

#define MIDDLE 0
#define DOWN -1
#define UP 1
#define LEFT -1
#define RIGHT 1

#define MATRIX(x, y, n) (x*n + y)
#define CENTEROFMASSX(x, y) (grid.centerOfMassX[x*params.ncside + y])
#define CENTEROFMASSY(x, y) (grid.centerOfMassY[x*params.ncside + y])
#define MASS(x, y) (grid.m[x*params.ncside + y])

typedef struct _vector2{
	float x;
	float y;
} vector2;

typedef struct _vector2grid{
	int x;
	int y;
} vector2grid;

typedef struct _particle_t {
	int active;
	float m;
	double positionX;
	double positionY;
	double vx;
	double vy;
	long gridCoordinateX;
	long gridCoordinateY;
	double appliedForceX;
	double appliedForceY;
} particle_t;

typedef struct _grid_t{
	double *m;
	double *centerOfMassX;
	double *centerOfMassY;
} grid_t;

typedef struct _grid_tt{
	double m;
	double centerOfMassX;
	double centerOfMassY;
} grid_tt;

#define vectorNorm(x, y) (sqrt(x * x + y * y))
#define constrain1(size, n) ((n) == (-1) ? ((n) + (size)) : (n))
#define constrain2(size, n) ((n) == (size) ? ((n) - (size)) : (n))
#define SUM_A( x, y )  ((x) == 0 || (y) == 0 ? 0 : ( ( ( (x) * (x) ) / ( ( x ) + ( y ) ) ) * ( y ) ))
#define compareVectorsGrid(ax, ay, bx, by) ((ax) == (bx) && (ay) == (by) ? 1 : 0)


//int constrain(int size, int n);
/*float addVectors(float a, float b);
vector2 subVectors(vector2 a, vector2 b);
vector2 multiplyVectorByConst(float c, vector2 v);

vector2 calculateNextPosition(particle_t particle); // x = x0 + v0t + 0.5 a t^2 (t = 1)

vector2 calculateNextVelocity(particle_t particle); // v = v0 + at (t = 1)*/
void calculateGravForce(particle_t *p1, double massCenterX, double massCenterY, float m, int sideUPDOWN, int sideLEFTRIGHT);

#endif