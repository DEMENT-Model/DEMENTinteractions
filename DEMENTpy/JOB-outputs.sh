#!/bin/bash

#SBATCH -A allisons_lab           ## account to charge
#SBATCH -p free   	          ## run on the standard partition
#SBATCH -N 1                      ## run on a single node
#SBATCH --ntasks-per-node=1	  ## number of tasks to launch per node (standard x 4.5, highmem x 6)
#SBATCH --error=slurm-%J.err	  ## write errors in slurm-<jobID>.err file

module purge
hostname
module load python/3.8.0
cd src/

# usage: $1 relative path to output folder (../output-20230607)
python createCSVs.py  $1

# Move summary outputs to a different folder
cd $1
mkdir summaries
mv *.csv summaries
