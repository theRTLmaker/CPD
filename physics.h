#ifndef PHYSICS_H
#define PHYSICS_H

#define G 6.67408e-11

typedef struct _vector2{
	double x;
	double y;
} vector2;

typedef struct _particle_t{
	vector2 position;
	vector2 velocity;
	/*double x;
	double y;
	double vx;
	double vy;*/
	double m;
} particle_t;



void init_particles(long seed, long ncside, long long n_part, particle_t *par);

vector2 addVectors(vector2 a, vector2 b);
vector2 subVectors(vector2 a, vector2 b);
vector2 multiplyVectorByConst(double c, vector2 v);

particle_t calculateCenterOfMass(particle_t *par, long long n_part);

vector2 calculateGravForce(particle_t p1, particle_t massCenter);

void printParticle(particle_t p);

#endif