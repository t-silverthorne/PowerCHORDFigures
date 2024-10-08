#!/bin/bash
#SBATCH --account=def-stinch
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000M               # memory per node
#SBATCH --time=0-00:10

nvidia-smi

# Choose a version of MATLAB by loading a module:
module load matlab/2023b.2
# Remove -singleCompThread below if you are using parallel commands:
matlab -nojvm -singleCompThread -batch "gpu_test"