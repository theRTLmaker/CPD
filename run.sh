#!/bin/bash
#SBATCH --job-name=--
#SBATCH --output=%j.txt
#SBATCH --error=%j.err
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=4
#SBATCH --nodes=32
srun ./simpar 5 200 300000000 2