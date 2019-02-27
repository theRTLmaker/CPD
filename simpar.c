#include "physics.h"
#include <stdlib.h>
#define RND0_1 ((double) random() / ((long long)1<<31))
#define EPSLON 0.01

void init_particles(long seed, long ncside, long long int n_part, particle_t *par)
{
    long long i;

    srandom(seed);

    for(i = 0; i < n_part; i++)
    {
        par[i].position.x = RND0_1;
        par[i].position.y = RND0_1;
        par[i].velocity.x = RND0_1 / ncside / 10.0;
        par[i].velocity.y = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}
