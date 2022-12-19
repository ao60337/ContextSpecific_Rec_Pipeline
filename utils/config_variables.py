from .medium_variables import UPTAKE_DRAINS
from numpy import log, linspace

MODEL = 'HumanGEM'
DATASET = 'Desai-GTEx'
# DATASET = 'Blanco-Melo'
ALGORITHMS = ['fastcore', 'tinit']
# ALGORITHMS = ['fastcore']
THRESHOLDS = linspace(0.25, 0.75, 5)
FASTCORE_THRESHOLDS = THRESHOLDS*2*log(2)
TINIT_THRESHOLDS = THRESHOLDS*2*log(2)


if DATASET == 'Desai-GTEx':
    THREAD_NUMBER_FASTCORE = 5
    THREAD_NUMBER_TINIT = 2

elif DATASET == 'Blanco-Melo':
    THREAD_NUMBER_FASTCORE = 4
    THREAD_NUMBER_TINIT = 2

PROTECTED = [['biomass_human', 'VBOF', 'EX_VBOF']] + list(UPTAKE_DRAINS['HAM Medium'])

