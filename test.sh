#!/bin/bash

time mpirun -np 4 ./simpar 1 3 10 2
time mpirun -np 4 ./simpar 3 10 100 5
time mpirun -np 4 ./simpar 3 5 500 5
time mpirun -np 4 ./simpar 2 5 10000 10
time mpirun -np 4 ./simpar 4 20 100000 5
time mpirun -np 4 ./simpar 8 5 1000000 50
time mpirun -np 4 ./simpar 2 20 100000000 5