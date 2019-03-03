#include <stdio.h>
#include "physics.h"
#include "debug.h"




void printGrid(grid_t grid, int size){
	for(int i = 0; i < size ; i++){
		for(int j = 0; j < size; j++){
			printf("Cell [%d][%d]: \n", i, j);
			printf("centro de massa: ");
			printParticle(grid.cells[i][j].massCenter);
			
			printParticleList(grid.cells[i][j].particles);
		}
	}
}

void printParticleList(particle_t *head){
	particle_t *cur = head;
	while(cur != NULL){
		printParticle(*cur);
		cur = cur-> nextParticle;
	}
}

