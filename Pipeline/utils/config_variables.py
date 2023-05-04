from .medium_variables import UPTAKE_DRAINS
# from numpy import linsp
# ace

# Model and dataset names ----------------------------------------------------------------------------------------------
MODEL = 'HumanGEM'
DATASET = 'Desai-GTEx'
OBJECTIVE = 'biomass_human'
MEDIUM_NAME = 'HAM Medium'
# ----------------------------------------------------------------------------------------------------------------------

# Omics Parameters -----------------------------------------------------------------------------------------------------
NOMENCLATURE = 'ensemble_gene_id'
OMICS_TYPE = 'transcriptomics'
THRESHOLDING_STRATEGY = 'local2'  # default, global, local1, local2
GLOBAL_THRESHOLD_LOWER = 1  # int
GLOBAL_THRESHOLD_UPPER = 3  # int
LOCAL_THRESHOLD = 2  # int
'''---Thresholding parameters---
Thresholding Strategies:
 - default: Does not use any thresholding strategy to filter the omics data;
 - global: Uses only a global threshold;
 - local1: Uses a global and a local threshold;
 - local2: Uses two global thresholds, one lower and one upper, and a local threshold.
The numbers in the thresholding options represent the position of the value to use. Currently the options are:
[0.1, 0.25, 0.5, 0.75, 0.9]; the threshold value will then be value of the dataset at that corresponding quantile.'''
# ----------------------------------------------------------------------------------------------------------------------

# Troppo parameters ----------------------------------------------------------------------------------------------------
ALGORITHMS = ['fastcore']
THREAD_NUMBER_FASTCORE = 10  # int.
THREAD_NUMBER_TINIT = 1  # int.
AND_OR_FUNCS = (min, sum)  # (min, max) or (min, sum).
INTEGRATION_THRESHOLDS = [0]  # List of ints or floats.
PROTECTED = ['biomass_human', 'VBOF', 'EX_VBOF'] + list(UPTAKE_DRAINS['HAM Medium'])  # List of Reactions to protect.
# ------------------------------------------------ ---------------------------------------------------------------------

# Task evaluation parameters -------------------------------------------------------------------------------------------
EVALUATE_TASKS = False
# ----------------------------------------------------------------------------------------------------------------------

# Gap-filling parameters -----------------------------------------------------------------------------------------------
GAP_FILLING = False
# ----------------------------------------------------------------------------------------------------------------------

# Model Reconstruction -------------------------------------------------------------------------------------------------
RECONSTRUCT_MODELS = True
# ----------------------------------------------------------------------------------------------------------------------
