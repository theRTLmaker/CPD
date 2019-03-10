#include <stdio.h>
#include <stdlib.h>
#include "physics.h"
#include "init_program.h"
#include "linkedList.h"
#include "debug.h"


/* Fazer um handler para tratar dos argumentos de entrada e passar para o init_particles*/

int main(int argc, char *argv[])
{
	FILE *fp;
	grid_t grid;
	vector2 centerOfMass;
	particle_t *par;
	parameters params; 
	int i = 0;
	int j = 0, m, k, aux1, aux2;

	fp = fopen("positions.txt", "w");

	par = handler_input(argc, argv, &params);
	
	// Init Grid
	grid = initGrid(grid, params.ncside);

	// Time Step simulation
	for (k = 0; k < params.timeStep; k++) {
		printf("TIME STEP %d\n", k);
		// Run throw all the cells and computes all the center of mass

		for (int i = 0; i < params.ncside; ++i) {
			for (int j = 0; j < params.ncside; ++j) {
				grid.m[i][j] = 0;
				grid.centerOfMass[i][j].x = 0;
				grid.centerOfMass[i][j].y = 0;
			}
		}

		// Calculate center of mass of each grid cell
		for(i = 0; i < params.n_part; i++) {
			grid.centerOfMass[par[i].gridCoordinate.x][par[i].gridCoordinate.y] = addVectors(multiplyVectorByConst(par[i].m, par[i].position) ,grid.centerOfMass[par[i].gridCoordinate.x][par[i].gridCoordinate.y]);
			grid.m[par[i].gridCoordinate.x][par[i].gridCoordinate.y] += par[i].m;
		}
		for (int i = 0; i < params.ncside; ++i){
			for (int j = 0; j < params.ncside; ++j){
				grid.centerOfMass[i][j] = multiplyVectorByConst(1/grid.m[i][j] ,grid.centerOfMass[i][j]);
			}
		}

		// Compute interactions
		// Run all particles
		for(i=0; i < params.n_part; i++){
			//printf("	PARTICLE %d\n", i);
			par[i].appliedForce.x = 0;
			par[i].appliedForce.y = 0;

			// Run the adjacent grids
			for(j = -1; j <= 1; j++) {
				for (m = -1; m <= 1; m++) {

					int sideUPDOWN = MIDDLE;
					int sideLEFTRIGHT = MIDDLE;

					if(par[i].gridCoordinate.x+j == -1) 
						sideLEFTRIGHT = LEFT;
					if(par[i].gridCoordinate.x+j == params.ncside)
						sideLEFTRIGHT = RIGHT;
					if(par[i].gridCoordinate.y+m == -1)
						sideUPDOWN = DOWN;
					if(par[i].gridCoordinate.y+m == params.ncside)
						sideUPDOWN = UP;

					aux1 = constrain(params.ncside, par[i].gridCoordinate.x+j);

					aux2 = constrain(params.ncside, par[i].gridCoordinate.y+m);

					par[i].appliedForce = addVectors(par[i].appliedForce, calculateGravForce(par[i], grid.centerOfMass[ aux1 ][ aux2 ], grid.m[ aux1 ][ aux2 ], sideUPDOWN, sideLEFTRIGHT)); //for each adjacent cell.---- 
				}
			}
			// Updates particles position and velocity and position on the grid
			par[i].velocity = calculateNextVelocity(par[i]);
			par[i].position = calculateNextPosition(par[i]);

			par[i].pastPositions[k] = par[i].position;

			par[i].gridCoordinate = findPosition(par[i], params.ncside);
		}

	}
	

	centerOfMass.x = 0;
	centerOfMass.y = 0;
	long double totalMass = 0;
	for(int i = 0; i < params.n_part; i++) {
		centerOfMass = addVectors(centerOfMass, multiplyVectorByConst(par[i].m, par[i].position));
		totalMass += par[i].m;
	}
	centerOfMass = multiplyVectorByConst(1/totalMass, centerOfMass);


	//printGrid(grid, params.ncside);

	printf("Cell 0:\n x=%Lf\n y=%Lf\n", par[0].position.x, par[0].position.y);
	printf("Center of mass:\n x=%Lf\n y=%Lf\n", centerOfMass.x, centerOfMass.y);
	
	/*FILE output for debugging*/
	/*for (i = 0; i < params.n_part; i++)
	{
		for ( j = 0; j < params.timeStep; j++)
		{
			fprintf(fp, "%Lf,%Lf\r\n", par[i].pastPositions[j].x, par[i].pastPositions[j].y);
		}
		printf("----\n");
		free(par[i].pastPositions);
	}
	fclose(fp);*/
	/*FILE output for debugging*/

	freeEverything(par, grid, params.ncside);
	return 0;
}

