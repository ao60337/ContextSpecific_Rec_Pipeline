from .medium_variables import UPTAKE_DRAINS
from numpy import log, linspace

# Model and dataset names ----------------------------------------------------------------------------------------------
MODEL = 'Dsalina'
DATASET = 'nacl_h2o2_sorb' #'light'
OBJECTIVE = 'e_Biomass__cytop'
MEDIUM_NAME = 'f2 medium'
# ----------------------------------------------------------------------------------------------------------------------

# Omics Parameters -----------------------------------------------------------------------------------------------------
NOMENCLATURE = 'NCBI'
OMICS_TYPE = 'transcriptomics'
THRESHOLDING_STRATEGY = 'Local2'  # Default, Global, Local1, Local2
GLOBAL_THRESHOLD_UPPER = 2  # int
GLOBAL_THRESHOLD_LOWER = 4  # int
LOCAL_THRESHOLD = 4  # int
'''---Thresholding parameters---
Thresholding Strategies:
 - Default: Does not use any thresholding strategy to filter the omics data;
 - Global: Uses only a global threshold;
 - Local1: Uses a global and a local threshold;
 - Local2: Uses two global thresholds, one lower and one upper, and a local threshold.
The numbers in the thresholding options represent the position of the value to use. Currently the options are:
[0.1, 0.25, 0.5, 0.75, 0.9]; the threshold value will then be value of the dataset at that corresponding quantile.'''
# ----------------------------------------------------------------------------------------------------------------------

# Troppo parameters ----------------------------------------------------------------------------------------------------
ALGORITHMS = ['fastcore']
THREAD_NUMBER_FASTCORE = 10  # int.
THREAD_NUMBER_TINIT = 10  # int.
AND_OR_FUNCS = (min, sum)  # (min, max) or (min, sum).
# THRESHOLDS = linspace(0.25, 0.75, 5)  # List of ints or floats.
INTEGRATION_THRESHOLDS = [2.8, 3.7]  # NOT NECESSARY: Just a transformation of the threshold.
PROTECTED = ['e_Biomass__cytop'] + list(UPTAKE_DRAINS['f2 medium'])  # List of Reactions to protect.
# ----------------------------------------------------------------------------------------------------------------------

# Task evaluation parameters -------------------------------------------------------------------------------------------
EVALUATE_TASKS = False
# ----------------------------------------------------------------------------------------------------------------------

# Gap-filling parameters -----------------------------------------------------------------------------------------------
GAP_FILLING = False
# ----------------------------------------------------------------------------------------------------------------------

# Model Reconstruction -------------------------------------------------------------------------------------------------
RECONSTRUCT_MODELS = True
# ----------------------------------------------------------------------------------------------------------------------
