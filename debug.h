#ifndef DEBUG_H
#define DEBUG_H

void printVectorPosition(vector2 p);
void printVectorVelocity(vector2 v);
void printVectorGrid(vector2grid g);

void printParticle(particle_t p);

void printAllParticles(particle_t *p, long long nr_part);

void printCenter(particle_t p);

#endif