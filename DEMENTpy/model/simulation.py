"""
This module, simulation.py, contains the DEMENT() function for conducting the simulation.

"""

import numpy as np
import pandas as pd
import pickle

from grid import Grid
from output import Output

# Define the function which takes two inputs: a runtime file and a data_init file from  
def DEMENT(runtime_i, data_init_i):

  # A few system constraints
  pulse      = int(runtime_i.loc['pulse',1])         # number of pulses
  cycle      = int(runtime_i.loc['end_time',1])      # number of time steps in each pulse
  interval   = int(runtime_i.loc['interval',1])      # interval of time step to record outputs
  mic_reinit = True    # indicate reinitialization of microbial community
  
  # Prepare for output by creating an instance of the Output class
  Output_init = Output(runtime_i, data_init_i)
  
  # Create an instance of the Grid class
  Ecosystem = Grid(runtime_i, data_init_i)
  
  # Run the model 
  for p in range(pulse):
      
      for i in range(p*cycle, (p+1)*cycle):
      
          # substrates degradation
          Ecosystem.degradation(p,i)
      
          # monomers uptake
          Ecosystem.uptake(p,i)
      
          # microbial metabolism
          Ecosystem.metabolism(i)
      
          # microbial death
          Ecosystem.mortality(i)
      
          # microbial reproduction and dispersal
          Ecosystem.reproduction(i)
      
          # output data using the "output" method in the Output class
          if i == 0:
              Output_init.output(Ecosystem,i)  # day 1
          elif i%interval==interval-1:
              Output_init.output(Ecosystem,i)  # interval
          
          # record the grid at time points divisible by 40
          if i%40==0:
              Output_init.outputGrid(Ecosystem, i)
          
          # if only 1 pusle, skip all following lines within this loop
          #if pulse == 1:
          #    continue
          
          # output microbial mass of every iteration using the "microbes_df" method in the Output class
          Output_init.microbes_abundance(Ecosystem,i)
          
          # re-initialize microbial community in each new pulse
          if i == (p+1)*cycle-1:
              Ecosystem.repopulation(Output_init,i,mic_reinit)
  return Output_init
