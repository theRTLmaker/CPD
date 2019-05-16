#!/bin/bash
#SBATCH --job-name=---
#SBATCH --output=%j.txt
#SBATCH --error=%j.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --nodes=8
srun ./simpar 8 1500 200000000 5