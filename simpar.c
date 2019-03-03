#include <stdio.h>
#include <stdlib.h>
#include "physics.h"
#include "init_program.h"
#include "LinkedList.h"
#include "debug.h"

/* Fazer um handler para tratar dos argumentos de entrada e passar para o init_particles*/

int main(int argc, char *argv[])
{
	FILE *fp;
	grid_t grid;
	particle_t *par;
	particle_t center;
	parameters params; 
	int i = 0;
	int j = 0, m, k, aux1, aux2;

	fp = fopen("positions.txt", "w");

	par = handler_input(argc, argv, &params);
	
	//Init Grid
	grid = initGrid(grid, params.ncside);
	//Assign particles to corresponding cell
	for(i = 0; i < params.n_part; i++){
		grid.cells[ par[i].gridCoordinate.x ][ par[i].gridCoordinate.y ].particles = addParticle(grid.cells[ par[i].gridCoordinate.x ][ par[i].gridCoordinate.y ].particles, &par[i]);
	}


	for (k = 0; k < params.timeStep; k++){
		for(i = 0; i < params.ncside; i++){
			for (j = 0; j < params.ncside; j++)
			{
				grid.cells[i][j].massCenter = calculateCenterOfMass(grid.cells[i][j].particles); //For each cell computes center of mass
			}
		}

		for(i=0; i < params.n_part; i++){
			for(j = -1; j <= 1; j++){
				for (m = -1; m <= 1; m++){
					aux1 = constrain(par[i].gridCoordinate.x+j);

					aux2 = constrain(par[i].gridCoordinate.y+m);



					par[i].appliedForce = addVectors(par[i].appliedForce, calculateGravForce(par[i], grid.cells[ aux1 ][ aux2 ].massCenter)); //for each adjacent cell.---- 
				}
			}
			/*par[i].velocity = calculateNextVelocity(par[i]);
			par[i].position = calculateNextPosition(par[i]);
			par[i] = findPosition(par[i], params.ncside); //position in grid*/

			//if particle changed cell
			//		remove particle form cell list
			//		add to new cell
		}
	}
	printGrid(grid, params.ncside);
	
/*	for (j = 0; j < params.timeStep; j++)
	{
		for(i = 0; i<params.n_part; i++){
			//for(int k = 0; k < 9; k++){
				center = calculateCenterOfMass(par, params.n_part); 
				//Iterar por cada célula adjacente e a própria célula
				//calcular forças em relação aos centros de massa e somá-las
				par[i].appliedForce = calculateGravForce(par[i], center);
			//}
			
			par[i].velocity = calculateNextVelocity(par[i]);
			par[i].position = calculateNextPosition(par[i]);
			par[i].pastPositions[j] = par[i].position;
		}

		

		printAllParticles(par, params.n_part);
		
		printCenter(center);
	}*/


	/*FILE output for debugging*/
	for (i = 0; i < params.n_part; i++)
	{
		for ( j = 0; j < params.timeStep; j++)
		{
			fprintf(fp, "%Lf,%Lf\r\n", par[i].pastPositions[j].x, par[i].pastPositions[j].y);
		}
		free(par[i].pastPositions);
	}
	fclose(fp);
	/*FILE output for debugging*/

	freeEverything(par, grid.cells, params.ncside);
	return 0;
}

