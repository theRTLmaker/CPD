#ifndef PHYSICS_H
#define PHYSICS_H

#define G 6.67408e-11
#define EPSLON 0.01

#define MIDDLE -1
#define DOWN 0
#define UP 1
#define LEFT 2
#define RIGHT 3


typedef struct _vector2{
	long double x;
	long double y;
} vector2;

typedef struct _vector2grid{
	int x;
	int y;
} vector2grid;

typedef struct _particle_t {
	long double m;
	vector2 position;
	vector2 velocity;
	vector2grid gridCoordinate;
	vector2 appliedForce;
	vector2 *pastPositions;
	struct _particle_t * nextParticle;
} particle_t;

typedef struct _grid_t{
	long double **m;
	vector2 **centerOfMass;
} grid_t;

#define constrain(size, n) ((n)==-1 ? (n + size):(((n)==size?(n - size):n)):n)
#define SUM_A( x, y )  ((x) == 0 || (y) == 0 ? 0 : ( ( ( (x) * (x) ) / ( ( x ) + ( y ) ) ) * ( y ) ))


//int constrain(int size, int n);
vector2 addVectors(vector2 a, vector2 b);
vector2 subVectors(vector2 a, vector2 b);
vector2 multiplyVectorByConst(double c, vector2 v);

int compareVectorsGrid(vector2grid a, vector2grid b);

particle_t calculateCenterOfMass(particle_t *head);

vector2 calculateGravForce(particle_t p1, vector2 massCenter, long double m, int sideUPDOWN, int sideLEFTRIGHT);

vector2 calculateNextPosition(particle_t particle); // x = x0 + v0t + 0.5 a t^2 (t = 1)

vector2 calculateNextVelocity(particle_t particle); // v = v0 + at (t = 1)


#endif