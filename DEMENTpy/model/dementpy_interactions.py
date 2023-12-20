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

# Grow a seed of random number generator
np.random.seed(int(microbe_seed))

# Initialize the 20 taxa community
microbe_initialization = initialize_microbe(runtime, grid_initialization)

# Combine grid and microbes into a single initialization
data_init = initialize_data(runtime, grid_initialization, microbe_initialization)

# Run simulations
output=DEMENT(runtime, data_init)
print('simulation complete')
export(output,'grid'+grid_seed+'_microbe'+microbe_seed+'_all')
