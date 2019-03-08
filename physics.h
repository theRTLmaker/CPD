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

typedef struct _gridcell{
	particle_t massCenter;
	particle_t *particles;
} gridCell;

typedef struct _grid_t{
	long  dimension;
	//particle_t ***cells;
	gridCell **cells;
} grid_t;



int constrain(int size, int n);
vector2 addVectors(vector2 a, vector2 b);
vector2 subVectors(vector2 a, vector2 b);
vector2 multiplyVectorByConst(double c, vector2 v);

int compareVectorsGrid(vector2grid a, vector2grid b);

particle_t calculateCenterOfMass(particle_t *head);

vector2 calculateGravForce(particle_t p1, particle_t massCenter, int sideUPDOWN, int sideLEFTRIGHT);

vector2 calculateNextPosition(particle_t particle); // x = x0 + v0t + 0.5 a t^2 (t = 1)

vector2 calculateNextVelocity(particle_t particle); // v = v0 + at (t = 1)

void printParticle(particle_t p);

void printAllParticles(particle_t *p, long long int nr_part);

void printCenter(particle_t p);

#endif