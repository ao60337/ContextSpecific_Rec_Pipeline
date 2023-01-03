from .medium_variables import UPTAKE_DRAINS
from numpy import log, linspace

# Model and dataset names ----------------------------------------------------------------------------------------------
MODEL = 'HumanGEM'
# DATASET = 'Desai-GTEx'
# DATASET = 'Blanco-Melo'
DATASET = 'CCLE_expression_full'
OBJECTIVE = 'biomass_human'
MEDIUM_NAME = 'HAM Medium'
# ----------------------------------------------------------------------------------------------------------------------

# Omics Parameters -----------------------------------------------------------------------------------------------------
NOMENCLATURE = 'ensemble_gene_id'
OMICS_TYPE = 'transcriptomics'
THRESHOLDING_STRATEGY = 'Local2'  # Default, Global, Local1, Local2
GLOBAL_THRESHOLD_UPPER = 0  # int
GLOBAL_THRESHOLD_LOWER = 3  # int
LOCAL_THRESHOLD = 1  # int
# ----------------------------------------------------------------------------------------------------------------------

# Troppo parameters ----------------------------------------------------------------------------------------------------
# ALGORITHMS = ['fastcore', 'tinit']
ALGORITHMS = ['fastcore']
THREAD_NUMBER_FASTCORE = 10  # int.
THREAD_NUMBER_TINIT = 4  # int.
AND_OR_FUNCS = (min, sum)  # (min, max) or (min, sum).
THRESHOLDS = linspace(0.25, 0.75, 5)  # List of ints or floats.
INTEGRATION_THRESHOLDS = [THRESHOLDS[0] * 2 * log(2)]  # NOT NECESSARY: Just a transformation of the thresholds.
PROTECTED = [['biomass_human', 'VBOF', 'EX_VBOF']] + list(UPTAKE_DRAINS['HAM Medium'])  # List of Reactions to protect.
# ----------------------------------------------------------------------------------------------------------------------

# Task evaluation parameters -------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# Gap-filling parameters -----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
