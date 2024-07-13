#!/bin/bash -l
#SBATCH --job-name=matlab_sweep_test
#SBATCH --array=1-3
#SBATCH --account=def-stinch 
#SBATCH --time=0-00:05
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      
#SBATCH --mem-per-cpu=2000             
#SBATCH --mail-user=turner.silverthorne@gmail.com
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2023b.2
# Remove -singleCompThread below if you are using parallel commands:
matlab -singleCompThread -batch "sweep_yalmip"
