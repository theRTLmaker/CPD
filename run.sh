#!/bin/bash
#SBATCH --job-name=group13
#SBATCH --output=%j.txt
#SBATCH --error=%j.err
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
srun ./simpar 2 20 100000000 5