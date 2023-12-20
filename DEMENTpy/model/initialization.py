"""
This module, with only one function "initialized_data(), initializes data related to 
substrate, monomer, enzyme, and microbe, as well as their distribution on the spatial grid,
preceding the actual decompostion-related computations.
"""
import pandas as pd
import numpy as np

from substrate import Substrate
from monomer   import Monomer
from enzyme    import Enzyme
from microbe   import Microbe
from utility   import expand

def initialize_grid(runtime_parameters):
    """
    Parameters:
        runtime_parameters: user-specified parameters setting up the system;
                            all other paras loaded by reading the parameters.csv
    Return:
        Data_Dictionary: a dictionary of all variables of the substrates, monomers, and enzymes
    """
    # Load all input files
    parameters      = pd.read_csv('parameters.csv',         header=None, index_col=0).astype('float32')   # parameters
    substrates_init = pd.read_csv('initial_substrates.csv', header=0,    index_col=0).astype('float32')   # initial substrates
    sub_mon_input   = pd.read_csv('sub_mon_inputs.csv',     header=0,    index_col=0).astype('float32')   # inputs of substrates and monomers
    Ea_input        = pd.read_csv("enzyme_ea.csv",          header=0,    index_col=0).astype('float32')   # enzyme activation energy
    climate         = pd.read_csv('climate.csv',            header=0,    index_col=0)                     # climate forcings

    # daily temperature and water potential
    daily_temp = climate['Temp'].astype('float32')  # temperaure series
    daily_psi  = climate['Psi'].astype('float32')   # water potential series

    #...an instance of Substrate class 
    Substrates = Substrate(runtime_parameters,parameters,substrates_init)
    #...substrate initial pool size
    substrates_initial_pool = Substrates.Substrates_start
    #...substrate input rate
    substrates_input_rate = Substrates.substrate_input(sub_mon_input)
    #...substrates-produced monomers
    substrates_produced_monomers = Substrates.substrate_produced_monomer()
    #...substrates degradation required enzymes
    substrates_req_enzyme = Substrates.substrate_degradation_enzyme()
    
    #...an instance of Monomer class
    Monomers = Monomer(runtime_parameters,parameters)
    #...monomers initial pool size
    monomers_initial_pool = Monomers.monomer_initialization(substrates_initial_pool)
    #...initial monomer ratios
    monomer_ratio_inital = Monomers.monomer_ratios(monomers_initial_pool)
    #...monomers input rate
    monomers_input_rate = Monomers.monomer_input_rate(sub_mon_input)
    #...monomers uptake required enzymes
    monomers_uptake_reqenzyme = Monomers.monomer_uptake_reqenzyme()
    
    #...an instance of Enzyme class
    Enzymes = Enzyme(runtime_parameters,parameters,substrates_initial_pool.index)
    #...enzyme initial pool size:0
    enzymes_initial_pool = Enzymes.enzyme_pool_initialization()
    #...enzyme attributes
    enzymes_attributes = Enzymes.enzyme_attributes()
    #...enzymes of substrate degradation Ea    
    enzymes_Ea = Enzymes.enzyme_Ea(Ea_input)
    #...monomers uptake enzyme Ea
    enzymes_uptake_Ea = Enzymes.enzyme_uptake_Ea()
    #...enzymes of substrate degradation Vmax
    enzymes_Vmax,enzymes_Vmax_T = Enzymes.enzyme_Vmax(substrates_req_enzyme)
    #...monomers uptake enzyme Vmax
    enzymes_uptake_Vmax= Enzymes.enzyme_uptake_Vmax(monomers_uptake_reqenzyme)
    #...enzymes of substrate degradation Km
    enzymes_Km = Enzymes.enzyme_Km(enzymes_Vmax)
    #...monomers uptake enzyme Km
    enzymes_uptake_Km = Enzymes.enzyme_uptake_Km(enzymes_uptake_Vmax)
    
    #...Store the data in a dictionary
    Data_Dictionary = {"substrates_initial_pool": substrates_initial_pool,
                       "substrates_input_rate":  substrates_input_rate,
                       "substrates_produced_monomers": substrates_produced_monomers,
                       "substrates_req_enzyme":  substrates_req_enzyme,
                       "monomers_initial_pool": monomers_initial_pool,
                       "monomer_ratio_inital": monomer_ratio_inital,
                       "monomers_input_rate": monomers_input_rate,
                       "monomers_uptake_reqenzyme": monomers_uptake_reqenzyme,
                       "enzymes_initial_pool": enzymes_initial_pool,
                       "enzymes_attributes": enzymes_attributes,
                       "enzymes_Ea": enzymes_Ea,
                       "enzymes_uptake_Ea": enzymes_uptake_Ea,
                       "enzymes_Vmax_T": enzymes_Vmax_T,
                       "enzymes_uptake_Vmax": enzymes_uptake_Vmax,
                       "enzymes_Km": enzymes_Km,
                       "enzymes_uptake_Km": enzymes_uptake_Km
    }
    return Data_Dictionary

def initialize_microbe(runtime_parameters, grid_init):
    """
    Parameters:
        runtime_parameters: user-specified parameters setting up the system;
                            all other paras loaded by reading the parameters.csv
    Return:
        Data_Dictionary: a dictionary of all variables of the microbes 
    """
    # Load all input files
    parameters      = pd.read_csv('parameters.csv',         header=None, index_col=0).astype('float32')   # parameters
    substrates_init = pd.read_csv('initial_substrates.csv', header=0,    index_col=0).astype('float32')   # initial substrates
    sub_mon_input   = pd.read_csv('sub_mon_inputs.csv',     header=0,    index_col=0).astype('float32')   # inputs of substrates and monomers
    Ea_input        = pd.read_csv("enzyme_ea.csv",          header=0,    index_col=0).astype('float32')   # enzyme activation energy
    climate         = pd.read_csv('climate.csv',            header=0,    index_col=0)                     # climate forcings

    # daily temperature and water potential
    daily_temp = climate['Temp'].astype('float32')  # temperaure series
    daily_psi  = climate['Psi'].astype('float32')   # water potential series
    
    #...an instance of Microbe class
    Microbes = Microbe(runtime_parameters,parameters)
    #...Microbial community initialization#...note microbial_community is a tuple
    microbial_community = Microbes.microbial_community_initialization()
    #...Microbial minimum ratios
    microbial_min_ratios = Microbes.minimum_cell_quota()
    #...Microbial enzyme genes
    microbial_enzyme_gene = Microbes.microbe_enzyme_gene()
    #...Microbial osmolyte genes
    microbial_osmolyte_gene = Microbes.microbe_osmolyte_gene()
    #...Microbial uptake genes
    microbial_uptake_gene = Microbes.microbe_uptake_gene(grid_init["substrates_req_enzyme"],microbial_enzyme_gene,grid_init["substrates_produced_monomers"])
    #...Microbial uptake cost
    microbial_uptake_cost = Microbes.microbe_uptake_cost(microbial_uptake_gene)
    #...Microbial enzyme production rate
    microbial_enzyme_prod_rate = Microbes.microbe_enzproduction_rate(microbial_enzyme_gene,grid_init["enzymes_attributes"])
    #...Microbial osmolyte productoin rate
    microbial_osmolyte_prod_rate = Microbes.microbe_osmoproduction_rate(microbial_osmolyte_gene)
    #...Microbial drought tolerance
    microbial_drought_tol = Microbes.microbe_drought_tol(microbial_osmolyte_prod_rate[2],microbial_osmolyte_prod_rate[3])
    #...Microbial mortality
    microbial_mortality = Microbes.microbe_mortality(microbial_community[2])
    
    #...Store the data in a dictionary
    Data_Dictionary = {"microbial_community":microbial_community,
                       "microbial_min_ratios":  microbial_min_ratios,
                       "microbial_enzyme_gene": microbial_enzyme_gene,
                       "microbial_osmolyte_gene": microbial_osmolyte_gene,
                       "microbial_uptake_gene": microbial_uptake_gene,
                       "microbial_uptake_cost": microbial_uptake_cost,
                       "microbial_enzyme_prod_rate": microbial_enzyme_prod_rate,
                       "microbial_osmolyte_prod_rate": microbial_osmolyte_prod_rate, 
                       "microbial_drought_tol": microbial_drought_tol,
                       "microbial_mortality": microbial_mortality
    }
    return Data_Dictionary
    
def initialize_data(runtime_parameters, grid_init, microbe_init):
    """
    Parameters:
        runtime_parameters: user-specified parameters setting up the system;
                            all other paras loaded by reading the parameters.csv
    Return:
        Data_Dictionary: a dictionary of all variables that feeds the grid.py module
    """
    
    # Load all input files
    parameters      = pd.read_csv('parameters.csv',         header=None, index_col=0).astype('float32')   # parameters
    substrates_init = pd.read_csv('initial_substrates.csv', header=0,    index_col=0).astype('float32')   # initial substrates
    sub_mon_input   = pd.read_csv('sub_mon_inputs.csv',     header=0,    index_col=0).astype('float32')   # inputs of substrates and monomers
    Ea_input        = pd.read_csv("enzyme_ea.csv",          header=0,    index_col=0).astype('float32')   # enzyme activation energy
    climate         = pd.read_csv('climate.csv',            header=0,    index_col=0)                     # climate forcings

    # daily temperature and water potential
    daily_temp = climate['Temp'].astype('float32')  # temperaure series
    daily_psi  = climate['Psi'].astype('float32')   # water potential series

    #...substrate initial pool size
    substrates_initial_pool = grid_init['substrates_initial_pool']
    #...substrate input rate
    substrates_input_rate = grid_init['substrates_input_rate']
    #...substrates-produced monomers
    substrates_produced_monomers = grid_init['substrates_produced_monomers']
    #...substrates degradation required enzymes
    substrates_req_enzyme = grid_init['substrates_req_enzyme']
    
    #...monomers initial pool size
    monomers_initial_pool = grid_init['monomers_initial_pool']
    #...initial monomer ratios
    monomer_ratio_inital = grid_init['monomer_ratio_inital']
    #...monomers input rate
    monomers_input_rate = grid_init['monomers_input_rate']
    #...monomers uptake required enzymes
    monomers_uptake_reqenzyme = grid_init['monomers_uptake_reqenzyme']
    
    #...enzyme initial pool size:0
    enzymes_initial_pool = grid_init['enzymes_initial_pool']
    #...enzyme attributes
    enzymes_attributes = grid_init['enzymes_attributes']
    #...enzymes of substrate degradation Ea    
    enzymes_Ea = grid_init['enzymes_Ea']
    #...monomers uptake enzyme Ea
    enzymes_uptake_Ea = grid_init['enzymes_uptake_Ea']
    #...enzymes of substrate degradation Vmax
    enzymes_Vmax_T = grid_init['enzymes_Vmax_T']
    #...monomers uptake enzyme Vmax
    enzymes_uptake_Vmax= grid_init['enzymes_uptake_Vmax']
    #...enzymes of substrate degradation Km
    enzymes_Km = grid_init['enzymes_Km']
    #...monomers uptake enzyme Km
    enzymes_uptake_Km = grid_init['enzymes_uptake_Km']
    
    #...Microbial community initialization#...note microbial_community is a tuple
    microbial_community = microbe_init['microbial_community']
    #...Microbial minimum ratios
    microbial_min_ratios = microbe_init['microbial_min_ratios']
    #...Microbial enzyme genes
    microbial_enzyme_gene = microbe_init['microbial_enzyme_gene']
    #...Microbial osmolyte genes
    microbial_osmolyte_gene = microbe_init['microbial_osmolyte_gene']
    #...Microbial uptake genes
    microbial_uptake_gene = microbe_init['microbial_uptake_gene']
    #...Microbial uptake cost
    microbial_uptake_cost = microbe_init['microbial_uptake_cost']
    #...Microbial enzyme production rate
    microbial_enzyme_prod_rate = microbe_init['microbial_enzyme_prod_rate']
    #...Microbial osmolyte productoin rate
    microbial_osmolyte_prod_rate = microbe_init['microbial_osmolyte_prod_rate']
    #...Microbial drought tolerance
    microbial_drought_tol = microbe_init['microbial_drought_tol']
    #...Microbial mortality
    microbial_mortality = microbe_init['microbial_mortality']
    
    #...Dump all initialized data into a dictionary; NOTE: variables with expand() put on the spatial grid
    gridsize = int(runtime_parameters.loc['gridsize',1])
    
    Data_Dictionary = {"Substrates": expand(substrates_initial_pool,gridsize),
                       "SubInput":   expand(substrates_input_rate,gridsize),
                       "ReqEnz":           substrates_req_enzyme,
                       "MonomersProduced": substrates_produced_monomers,
                       "Monomers":     expand(monomers_initial_pool,gridsize),
                       "Monomer_ratio":expand(monomer_ratio_inital,gridsize),
                       "MonInput":     expand(monomers_input_rate,gridsize),
                       "Uptake_ReqEnz":expand(monomers_uptake_reqenzyme,gridsize),
                       "Enzymes":      expand(enzymes_initial_pool,gridsize),
                       "Km0":          expand(enzymes_Km,gridsize),            # enzyme half-saturation constant
                       "Uptake_Km0":   expand(enzymes_uptake_Km,gridsize),     # transporter half-saturation constant
                       "Uptake_Ea":    expand(enzymes_uptake_Ea,gridsize),     # transporter acitivation energy
                       "Uptake_Vmax0": expand(enzymes_uptake_Vmax,gridsize),   # transporter reaction rate
                       "Ea":           expand(enzymes_Ea,gridsize),            # enzyme activation energy
                       "Vmax0":        expand(enzymes_Vmax_T,gridsize),        # enzyme reaction rate
                       "EnzAttrib":    enzymes_attributes,                     # enzyme stoichiometry and energy cost
                       "Microbes_pp": microbial_community[0],                  # tuple[0]: microbes preceding placement
                       "Microbes":    microbial_community[1],                  # tuple[1]: initialized spatial microbes
                       "fb":          microbial_community[2],                  # tuple[2]: fungi index
                       "Bac_density": microbial_community[3],                  # tuple[3]: bacterial density
                       "Fun_density": microbial_community[4],                  # tuple[4]: fungi density
                       "MinRatios":   expand(microbial_min_ratios,gridsize),   # microbial cell min. ratios
                       "UptakeGenes": expand(microbial_uptake_gene,gridsize),  # transporter gene distribution across taxa
                       "OsmoGenes":   expand(microbial_osmolyte_gene,gridsize),# osmolyte gene distribution across taxa
                       "EnzGenes":    expand(microbial_enzyme_gene,gridsize),  # enzyme gene distribution across taxa
                       "UptakeGenes_trait":   expand(microbial_uptake_cost[0],gridsize),        # single gene cost of transporter
                       "OsmoProdConsti_trait":expand(microbial_osmolyte_prod_rate[0],gridsize), # single gene cost of constitutive osmolyte
                       "OsmoProdInduci_trait":expand(microbial_osmolyte_prod_rate[1],gridsize), # single gene cost of inducible osmolyte
                       "EnzProdConsti_trait": expand(microbial_enzyme_prod_rate[0],gridsize),   # single gene cost of constitutive enzyme
                       "EnzProdInduci_trait": expand(microbial_enzyme_prod_rate[1],gridsize),   # single gene cost of inducible enzyme
                       "UptakeGenesCost":     expand(microbial_uptake_cost[1],gridsize),        # distribution of transporter gene cost across taxa
                       "OsmoProdConsti":      expand(microbial_osmolyte_prod_rate[2],gridsize), # distribution of consti. osmolyte gene cost across taxa
                       "OsmoProdInduci":      expand(microbial_osmolyte_prod_rate[3],gridsize), # distribution of induci. osmolyte gene cost across taxa
                       "EnzProdConstit":      expand(microbial_enzyme_prod_rate[2],gridsize),   # distribution of consti. enzyme gene cost across taxa
                       "EnzProdInduce":       expand(microbial_enzyme_prod_rate[3],gridsize),   # distribution of induci. enzyme gene cost across taxa
                       "TaxDroughtTol":       expand(microbial_drought_tol,gridsize),           # distribution of taxon-specific drought tol.
                       'basal_death_prob':  microbial_mortality[0],                # basal death probability
                       'death_rate':        microbial_mortality[1],                # change rate of death prob. agaist mositure
                       "AE_ref":            parameters.loc["CUE_ref",1],           # Reference assimilation efficiency: 0.5
                       "AE_temp":           parameters.loc["CUE_temp",1],          # AE temperature sensitivity; default: -0.016
                       'Uptake_Maint_cost': parameters.loc['Uptake_Maint_cost',1], # constant of transporter maintenence cost
                       'C_min':             parameters.loc['C_min',1],             # C threshold of cell lysis
                       'N_min':             parameters.loc['N_min',1],             # N threshold of cell lysis
                       'P_min':             parameters.loc['P_min',1],             # P threshold of cell lysis
                       'max_size_b':        parameters.loc['max_size_b',1],        # C quota threshold for bacterial cell division
                       'max_size_f':        parameters.loc['max_size_f',1],        # C quota threshold for fungal cell division
                       'wp_fc':             parameters.loc['wp_fc',1],             # threshold below which microbes start to respond to drought
                       'wp_th':             parameters.loc['wp_th',1],             # threshold below which microbes in full swing to respond to drought
                       'alpha':             parameters.loc['alpha',1],             # factor delineating curve concavity of microbial response to drought
                       'Temp': daily_temp,                                         # temperature
                       'Psi':  daily_psi                                           # water potential
                      }

    return Data_Dictionary
    
# Define function for sampling an initialization from the previous
def sample_microbe_init(microbe_initialization, taxa, n_taxa, taxa_per_box):

  # A few system constants
  runtime    = pd.read_csv('runtime.txt',header=None,index_col=0,sep='\t')
  gridsize = int(runtime.loc['gridsize',1])
  # Set up some rules
  # taxa=[0,2,3]
  # n_taxa=3
  # taxa_per_box=0.33
  max_size_b=2
  max_size_f=50
  fb = np.random.choice([1,0], n_taxa, replace=True, p=[runtime.loc['fb',1],(1-runtime.loc['fb',1])]).astype('int8') #index of fungal taxa in a microbial pool (1);1-d array
  BacC = 0.5  * max_size_b
  FunC = 0.5 * max_size_f


  ### Microbial Community List 
  # Microbes_pp
  df=microbe_initialization['microbial_community'][0].iloc[taxa,:]
  microbes_array=np.tile(df,(gridsize,1))
  index = ["Tax" + str(i) for i in [x+1 for x in taxa]] * gridsize
  microbes_pp = pd.DataFrame(data=microbes_array, index=index, columns=["C","N","P"], dtype='float32')

  # Export the microbial community preceding placement on the grid
  microbes_df = microbes_pp.copy(deep=True)

  # Derive the Fungi index by expanding the fb to the spatial grid
  fb_grid = np.tile(fb,gridsize)
  # Randomly place the microbial community created above on the spatial grid to initialize a spatially explicit microbial community by
  pb = taxa_per_box
  choose_taxa = np.random.choice([1,0], n_taxa*gridsize,replace=True, p=[pb,(1-pb)])
  pf = pb * max_size_b/max_size_f
  choose_taxa[fb_grid==1] = np.random.choice([1,0], sum(fb_grid),replace=True, p=[pf,(1-pf)])
  microbes_df.loc[choose_taxa==0] = 0

  # Derive the number of fungi and bacteria taxa
  Bac_index = microbes_df['C'] == BacC
  Fun_index = microbes_df['C'] == FunC
  Bac_taxa = microbes_df[Bac_index].groupby(level=0).sum().shape[0]
  Fun_taxa = microbes_df[Fun_index].groupby(level=0).sum().shape[0]
  print('After placement--','Bac_taxa=',Bac_taxa,'Fun_taxa=',Fun_taxa)

  Bac_density = sum(Bac_index)*BacC/gridsize
  Fun_density = sum(Fun_index)*FunC/gridsize
  print('After placement--','Bac_density=',Bac_density,'Fun_density=',Fun_density)

  # Place together
  microbial_community=[microbes_pp, microbes_df, fb_grid,Bac_density,Fun_density]

  ### Binary dataframes  
  microbial_min_ratios=microbe_initialization['microbial_min_ratios'].iloc[taxa,:]
  microbial_enzyme_gene=microbe_initialization['microbial_enzyme_gene'].iloc[taxa,:]
  microbial_osmolyte_gene=microbe_initialization['microbial_osmolyte_gene'].iloc[taxa,:] # 4 X 20 osmolytes binary dataframe
  microbial_uptake_gene=microbe_initialization['microbial_uptake_gene'].iloc[taxa,:] # 4 by 14 monomer binary dataframe

  ### Microbial_uptake_cost 
  microbial_uptake_cost=[microbe_initialization['microbial_uptake_cost'][0][taxa], microbe_initialization['microbial_uptake_cost'][1].iloc[taxa,:]]

  ### Microbial_enzyme_prod_rate
  microbial_enzyme_prod_rate=[microbe_initialization['microbial_enzyme_prod_rate'][0][taxa], microbe_initialization['microbial_enzyme_prod_rate'][1][taxa], microbe_initialization['microbial_enzyme_prod_rate'][2].iloc[taxa,:], microbe_initialization['microbial_enzyme_prod_rate'][3].iloc[taxa,:]]

  ### Microbial_osmolyte_prod_rate
  microbial_osmolyte_prod_rate=[microbe_initialization['microbial_osmolyte_prod_rate'][0][taxa], microbe_initialization['microbial_osmolyte_prod_rate'][1][taxa], microbe_initialization['microbial_osmolyte_prod_rate'][2].iloc[taxa,:], microbe_initialization['microbial_osmolyte_prod_rate'][3].iloc[taxa,:]]

  ### Microbial_drought_tol
  microbial_drought_tol=microbe_initialization['microbial_drought_tol'][taxa]

  ### Microbial_mortality
  microbial_mortality=[microbe_initialization['microbial_mortality'][0][0:n_taxa*gridsize], microbe_initialization['microbial_mortality'][1][0:n_taxa*gridsize]]

  ### Set up data dictionary
  Data_Dictionary = {"microbial_community":microbial_community,
                    "microbial_min_ratios":  microbial_min_ratios,
                    "microbial_enzyme_gene": microbial_enzyme_gene,
                    "microbial_osmolyte_gene": microbial_osmolyte_gene,
                    "microbial_uptake_gene": microbial_uptake_gene,
                    "microbial_uptake_cost": microbial_uptake_cost,
                    "microbial_enzyme_prod_rate": microbial_enzyme_prod_rate,
                    "microbial_osmolyte_prod_rate": microbial_osmolyte_prod_rate, 
                    "microbial_drought_tol": microbial_drought_tol,
                    "microbial_mortality": microbial_mortality
  }

  return Data_Dictionary

