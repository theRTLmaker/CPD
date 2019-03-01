#ifndef PHYSICS_H
#define PHYSICS_H

#define G 6.67408e-11

typedef struct _vector2{
	double x;
	double y;
} vector2;

typedef struct _vector2grid{
	int x;
	int y;
} vector2grid;

typedef struct _particle_t {
	double m;
	vector2 position;
	vector2 velocity;
	vector2grid gridCoordinate;
	//_particle_t *nextParticle;
} particle_t;

typedef struct _grid_t{
	long  dimension;
	particle_t **cells;
} grid_t;




vector2 addVectors(vector2 a, vector2 b);
vector2 subVectors(vector2 a, vector2 b);
vector2 multiplyVectorByConst(double c, vector2 v);

particle_t calculateCenterOfMass(particle_t *par, long long n_part);

vector2 calculateGravForce(particle_t p1, particle_t massCenter);

void printParticle(particle_t *p, long long int nr_part);

void printAllParticles(particle_t *p, long long int nr_part);

void printCenter(particle_t p);

#endif