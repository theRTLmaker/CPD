#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "physics.h"
#include "init_program.h"
#include "linkedList.h"
#include "debug.h"
#include <omp.h>
#include <time.h>


/* Fazer um handler para tratar dos argumentos de entrada e passar para o init_particles*/

int main(int argc, char *argv[])
{ 
	grid_t grid;
	particle_t *par;
	parameters params; 
	int k = 0;
	int sideUPDOWN;
	int sideLEFTRIGHT;
	int aux1 = 0, aux2 = 0;
	double invM;
	// init particles
	par = handler_input(argc, argv, &params);

    
	// Init Grid
	grid = initGrid(grid, params.ncside);


	// Time Step simulation
	for (k = 0; k < params.timeStep; k++) {
		
		// Run throw all the cells and resets all the center of mass
		memset(grid.m, 0, params.ncside*params.ncside*sizeof(double));
		memset(grid.centerOfMassX, 0, params.ncside*params.ncside*sizeof(double));
		memset(grid.centerOfMassY, 0, params.ncside*params.ncside*sizeof(double));
		
		
		// Calculate center of mass of each grid cell
		for(int i = params.n_part - 1; i >= 0; i--) {
			grid.centerOfMassX[MATRIX(par[i].gridCoordinateX, par[i].gridCoordinateY, params.ncside)] = grid.centerOfMassX[MATRIX(par[i].gridCoordinateX, par[i].gridCoordinateY, params.ncside)] + par[i].m * par[i].positionX;
			grid.centerOfMassY[MATRIX(par[i].gridCoordinateX, par[i].gridCoordinateY, params.ncside)] = grid.centerOfMassY[MATRIX(par[i].gridCoordinateX, par[i].gridCoordinateY, params.ncside)] + par[i].m * par[i].positionY;
			grid.m[MATRIX(par[i].gridCoordinateX, par[i].gridCoordinateY, params.ncside)] = grid.m[MATRIX(par[i].gridCoordinateX, par[i].gridCoordinateY, params.ncside)] + par[i].m;
		}

		for (int i = params.ncside*params.ncside - 1; i >= 0; i--)
		{
			grid.centerOfMassX[i] = grid.centerOfMassX[i]/grid.m[i];
			grid.centerOfMassY[i] = grid.centerOfMassY[i]/grid.m[i];
		}


		
		// Compute interactions
		// Run all particles
		for(int i = params.n_part - 1; i >= 0; i = i - 1){
			//printf("tid %d - PARTICLE %d\n", tid, i);
			par[i].appliedForceX = 0;
			par[i].appliedForceY = 0;

			// Run the adjacent grids
			int j;
			int m;
			for (int n = 0; n < 9; n = n + 1) {
				j = n / 3 - 1;
				m = n % 3 - 1;
				aux1 = 0, aux2 = 0;

				if(par[i].gridCoordinateX+j == -1) {
					sideLEFTRIGHT = LEFT;
					aux1 = par[i].gridCoordinateX+j + params.ncside;
				} else if(par[i].gridCoordinateX+j == params.ncside) {
					sideLEFTRIGHT = RIGHT;
					aux1 = par[i].gridCoordinateX+j - params.ncside;
				} else {
					sideLEFTRIGHT = MIDDLE;
					aux1 = par[i].gridCoordinateX+j;
				}

				if(par[i].gridCoordinateY+m == -1){
					sideUPDOWN = DOWN;
					aux2 = par[i].gridCoordinateY+m + params.ncside;
				} else if(par[i].gridCoordinateY+m == params.ncside){
					sideUPDOWN = UP;
					aux2 = par[i].gridCoordinateY+m - params.ncside;
				} else {
					sideUPDOWN = MIDDLE;
					aux2 = par[i].gridCoordinateY+m;
				}
				
				
				if(grid.m[ MATRIX(aux1, aux2, params.ncside) ] != 0)
					calculateGravForce(&(par[i]), grid.centerOfMassX[MATRIX(aux1, aux2, params.ncside) ], 
											grid.centerOfMassY[ MATRIX(aux1, aux2, params.ncside) ], 
											grid.m[ MATRIX(aux1, aux2, params.ncside) ], sideUPDOWN, sideLEFTRIGHT); //for each adjacent cell.---- 
			}

			/*for(int j = -1; j <= 1; j++) {
				for (int m = -1; m <= 1; m++) {
					aux1 = 0, aux2 = 0;

					if(par[i].gridCoordinateX+j == -1) {
						sideLEFTRIGHT = LEFT;
						aux1 = par[i].gridCoordinateX+j + params.ncside;
					} else if(par[i].gridCoordinateX+j == params.ncside) {
						sideLEFTRIGHT = RIGHT;
						aux1 = par[i].gridCoordinateX+j - params.ncside;
					} else {
						sideLEFTRIGHT = MIDDLE;
						aux1 = par[i].gridCoordinateX+j;
					}

					if(par[i].gridCoordinateY+m == -1){
						sideUPDOWN = DOWN;
						aux2 = par[i].gridCoordinateY+m + params.ncside;
					} else if(par[i].gridCoordinateY+m == params.ncside){
						sideUPDOWN = UP;
						aux2 = par[i].gridCoordinateY+m - params.ncside;
					} else {
						sideUPDOWN = MIDDLE;
						aux2 = par[i].gridCoordinateY+m;
					}
					
					
					if(grid.m[ MATRIX(aux1, aux2, params.ncside) ] != 0)
						calculateGravForce(&(par[i]), grid.centerOfMassX[MATRIX(aux1, aux2, params.ncside) ], 
												grid.centerOfMassY[ MATRIX(aux1, aux2, params.ncside) ], 
												grid.m[ MATRIX(aux1, aux2, params.ncside) ], sideUPDOWN, sideLEFTRIGHT); //for each adjacent cell.---- 
					
				}
			}*/
			invM = 1.0/par[i].m;
			// Updates particles position and velocity and position on the grid
			par[i].vx += par[i].appliedForceX * invM; //a = F/m 
			par[i].vy += par[i].appliedForceY * invM;

			par[i].positionX += par[i].vx + 0.5 * par[i].appliedForceX * invM;//x = x0 + v0t + 0.5 a t^2 (t = 1)
			par[i].positionY += par[i].vy + 0.5 * par[i].appliedForceY * invM;

			//See if its out of bounds
			if(par[i].positionX >= 1) par[i].positionX = par[i].positionX - (int)(par[i].positionX);
			else if(par[i].positionX < 0) par[i].positionX = 1 + (par[i].positionX - ceil(par[i].positionX)); 

			if(par[i].positionY >= 1) par[i].positionY = par[i].positionY - (int)(par[i].positionY);
			else if(par[i].positionY < 0) par[i].positionY = 1 + (par[i].positionY - ceil(par[i].positionY));

			// Updates the position of the particle on the grid of cells
			par[i].gridCoordinateX = par[i].positionX * params.ncside;
			par[i].gridCoordinateY = par[i].positionY * params.ncside;
			
		}
	}
	
	double centerOfMassX = 0;
	double centerOfMassY = 0;
	double totalMass = 0;

	// Calculates the center of mass of all cells
	for(int i = 0; i < params.n_part; i++) {
		centerOfMassX = centerOfMassX + par[i].m * par[i].positionX;
		centerOfMassY = centerOfMassY + par[i].m * par[i].positionY;
		totalMass = totalMass + par[i].m;
	}

	// prints the information
	centerOfMassX = centerOfMassX / totalMass;
	centerOfMassY = centerOfMassY / totalMass;
	printf("%.2f %.2f\n", par[0].positionX, par[0].positionY);
	printf("%.2f %.2f\n", centerOfMassX, centerOfMassY);

	freeEverything(par, grid, params.ncside);
	return 0;
}

