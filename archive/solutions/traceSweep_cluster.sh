#!/bin/bash
#SBATCH --account=def-stinch   # replace this with your own account
#SBATCH --cpus-per-task=32     # number of processes
#SBATCH --mem-per-cpu=750M     # memory; default unit is megabytes
#SBATCH --time=0-02:30         # time (DD-HH:MM)
#SBATCH --mail-user=turner.silverthorne@utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023
module load gcc/12.3
module load r/4.3.1  
module load gurobi/11.0.0

source ~/.profile
cd ~/research/powerCHORD2
Rscript solutions/traceSweep_cluster.R 