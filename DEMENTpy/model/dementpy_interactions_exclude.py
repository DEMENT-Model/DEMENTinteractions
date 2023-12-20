"""
-------------------------------------------------------------------------------
      DEMENTpy--Decomposition Model of Enzymatic Traits in Python,v1.0
                              Bin Wang, Ph.D.
              Department of Ecology and Evolutionary Biology
                       University of California Irvine
                  Emails: wbwenwu@gmail.com or bwang7@uci.edu
                          Twitter: @bio_atmosphere
-------------------------------------------------------------------------------
"""
import os
import sys
import pandas as pd
import numpy  as np
import pickle

from initialization import *
from grid import Grid
from output import Output
from utility import *
from simulation import DEMENT


print("""
  ---------------------------------------------------------------------------
         DEMENTpy (DEcomposition Model of Enzymatic Traits in Python)
                              Version 1.0
               Department of Ecology and Evolutionary Biology
                    University of California Irvine
  ---------------------------------------------------------------------------
  """)

# When trouble shooting:
# input_folder  = "input"  # input folder name
# output_folder = "output"   # output folder name
# outname       = "202206090930"  # output file name

# Obtain the command line arguments
input_folder  = sys.argv[1]   # input folder name
grid_seed     = sys.argv[2]   # output file name
microbe_seed  = sys.argv[3]
exclude       = sys.argv[4]   # taxa number ID to exclude
place_seed    = sys.argv[5]

# Set up the working directory
os.chdir('../'+input_folder)

# Grow a seed of random number generator
np.random.seed(int(grid_seed))

# A few system constants
runtime    = pd.read_csv('runtime.txt',header=None,index_col=0,sep='\t')
pulse      = int(runtime.loc['pulse',1])         # number of pulses
cycle      = int(runtime.loc['end_time',1])      # number of time steps in each pulse
interval   = int(runtime.loc['interval',1])      # interval of time step to record outputs
mic_reinit = True    # indicate reinitialization of microbial community

# Initialize data by calling the functions from initialize
grid_initialization = initialize_grid(runtime)

# Grow a seed of random number generator to initialize a community
np.random.seed(int(microbe_seed))

# Initialize the 25 taxa community
microbe_initialization = initialize_microbe(runtime, grid_initialization)

# Exclude the focal species
taxa_24=list(range(0, 25)) # create list of taxa IDs 0 through 25
taxa_24.remove(int(exclude)) # exclude a given taxa
microbe_init_24=sample_microbe_init(microbe_initialization, taxa=taxa_24, n_taxa=24,taxa_per_box=0.04) # create new initialization

# Build the runtime for only 19 taxa
runtime_24=build_runtime(n_taxa=24, taxa_per_box=0.04)

# Grow a seed of random number generator to change the placement of microbes on the grid
np.random.seed(int(place_seed))

# Combine grid and microbes into a single initialization
data_init_24 = initialize_data(runtime_24, grid_initialization, microbe_init_24)

# Run simulations
output_24=DEMENT(runtime_24, data_init_24)
print('simulation complete')
export(output_24,'grid'+grid_seed+'_microbe'+microbe_seed+'_exclude_'+exclude)
