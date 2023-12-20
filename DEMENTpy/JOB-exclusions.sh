#!/bin/bash

#SBATCH -A allisons_lab           ## account to charge
#SBATCH -p standard               ## run on the standard partition
#SBATCH -N 1                      ## run on a single node
#SBATCH --ntasks-per-node=1	  ## number of tasks to launch per node (standard x 4.5, highmem x 6)
#SBATCH --error=slurm-%J.err	  ## write errors in slurm-<jobID>.err file
#SBATCH --time=12:00:00           ## run time 12 hours

module purge
hostname
module load python/3.8.0
cd src

# python dementpy_interactions_exclude.py input grid_seed microbe_seed exclude place_seed

python dementpy_interactions_exclude.py $i 46819 $j $taxa $taxa$j

