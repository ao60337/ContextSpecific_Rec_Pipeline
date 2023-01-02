from medium_variables import UPTAKE_DRAINS
from numpy import log, linspace

# Model and dataset names ----------------------------------------------------------------------------------------------
MODEL = 'HumanGEM'
DATASET = 'Desai-GTEx_kidney_covid'
OBJECTIVE = 'biomass_human'
MEDIUM_NAME = 'HAM Medium'
# ----------------------------------------------------------------------------------------------------------------------

# Omics Parameters -----------------------------------------------------------------------------------------------------
NOMENCLATURE = 'ensemble_gene_id'
OMICS_TYPE = 'transcriptomics'
THRESHOLDING_STRATEGY = 'global'  # global, local1, local2
GLOBAL_THRESHOLD = 0  # int
LOCAL_THRESHOLD1 = 0  # int
LOCAL_THRESHOLD2 = 0  # int
# ----------------------------------------------------------------------------------------------------------------------

# Troppo parameters ----------------------------------------------------------------------------------------------------
ALGORITHMS = ['fastcore', 'tinit']
THREAD_NUMBER_FASTCORE = 10  # int.
THREAD_NUMBER_TINIT = 4  # int.
AND_OR_FUNCS = (min, max)  # (min, max) or (min, sum).
THRESHOLDS = linspace(0.25, 0.75, 5)  # List of ints or floats.
INTEGRATION_THRESHOLDS = THRESHOLDS * 2 * log(2)  # NOT NECESSARY: Just a transformation of the thresholds.
PROTECTED = [['biomass_human', 'VBOF', 'EX_VBOF']] + list(UPTAKE_DRAINS['HAM Medium'])  # List of Reactions to protect.
# ----------------------------------------------------------------------------------------------------------------------

# Task evaluation parameters -------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# Gap-filling parameters -----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
