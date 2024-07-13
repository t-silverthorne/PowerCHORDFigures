#!/bin/bash -l
#SBATCH --job-name=matlab_sweep_test
#SBATCH --array=1-228
#SBATCH --account=def-stinch 
#SBATCH --time=0-01:30
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-user=turner.silverthorne@gmail.com
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2023b.2
# Remove -singleCompThread below if you are using parallel commands:
matlab -nojvm -singleCompThread -batch "sweep_yalmip"
