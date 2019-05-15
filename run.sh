#!/bin/bash
#SBATCH --job-name=OMaisRapido
#SBATCH --output=%j.txt
#SBATCH --error=%j.err
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --nodes=16
srun ./simpar 8 1500 200000000 5