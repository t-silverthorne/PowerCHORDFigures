#!/bin/bash -l
#SBATCH --job-name=diff_ev_sweep
#SBATCH --array=1-280
#SBATCH --account=def-stinch 
#SBATCH --time=0-01:15
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000M               # memory per node
#SBATCH --mail-user=turner.silverthorne@gmail.com
#SBATCH --mail-type=ALL

nvidia-smi
# Choose a version of MATLAB by loading a module:
module load matlab/2023b.2
# Remove -singleCompThread below if you are using parallel commands:
matlab -nojvm -singleCompThread -batch "sweep_diffevolve"
