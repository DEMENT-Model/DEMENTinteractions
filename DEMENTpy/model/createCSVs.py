# Import modules
import os
import sys
import pandas as pd
import numpy  as np
import pickle

# Load custom modules
from output import Output

# Set working directory
os.chdir(sys.argv[1])

# Make empty respiration dataframe
carbonDF=pd.DataFrame()
averages=pd.DataFrame()
traitOutput=pd.DataFrame()
microbeDF=pd.DataFrame()
mortDF=pd.DataFrame()

# Use a for loop to open each pickle and extract different outputs
for files in os.listdir():
  if files.endswith(".pickle"):
    # Open a single pickle
    with open(files, 'rb') as f:
      x = pickle.load(f)

    #### Get timeseries of carbon contents
    ###############################################################
    # Make initial dataframe
    carboni=pd.DataFrame(list(range(0,1096)), columns = ['Time'])
    # Add file name
    carboni['file']=files
    # Add respiration
    carboni['respiration']=x.RespSeries
    # Get total biomass
    carboni['biomass']=x.MicrobesSeries.transpose().sum(axis=1)
    # Get total substrate
    carboni['substrate']=x.Substrates_Sum
    # Get total enzymes
    carboni['enzymes']=x.EnzymesSeries.transpose().sum(axis=1)
    # Append to the dataframe
    carbonDF=pd.concat([carbonDF, carboni])

    #### Get averages of carbon variables
    ###############################################################
    averagei={'file': files,
          'avgResp': carboni['respiration'].mean(),
          'avgBiomass':carboni['biomass'].mean(),
          'avgSubstrate':carboni['substrate'].mean()}
    averagei=pd.DataFrame(averagei, index=[0])
    # Append to the dataframe
    averages=pd.concat([averages, averagei])

    #### Get microbial time series
    ###############################################################
    microbei=x.MicrobesSeries.transpose()
    # Add file name
    microbei['file']=files
    # Add time column
    # microbei['Time']=microbei.index()
    # Append to the dataframe
    microbeDF=pd.concat([microbeDF, microbei])


    #### Get microbial time series
    ###############################################################
    morti=pd.DataFrame(list(range(0,1096)), columns = ['time'])
    # Add file name
    morti['file']=files
    # Add  death toll
    morti['biomassLoss']=x.biomassLoss
    # Append to the dataframe
    mortDF=pd.concat([mortDF, morti])


    #### Get trait tables
    ###############################################################
    if files.endswith("all.pickle"):
        # Create traitsi dataframe
    	traitsi={'file': files}
        traitsi=pd.DataFrame(traitsi, index=[0])
    	# Concat to microbial traits from Bin Wang
    	traitsi=pd.concat([traitsi, x.Microbial_traits], axis=1)
    	# Append to the dataframe
    	traitOutput=pd.concat([traitOutput, traitsi])

# Write out the csv
carbonDF.to_csv('carbonDF.csv', header=True)
averages.to_csv('averages.csv', header=True)
traitOutput.to_csv('traitOutput.csv', header=True)
microbeDF.to_csv('microbeDF.csv', header=True)
mortDF.to_csv('mortalityDF.csv', header=True)

