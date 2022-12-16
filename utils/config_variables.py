from .medium_variables import UPTAKE_DRAINS
from numpy import log, linspace

MODEL = 'HumanGEM'
DATASET = 'Desai-GTEx'
# DATASET = 'Blanco-Melo'
# ALGORITHMS = ['fastcore', 'tinit']
ALGORITHMS = ['tinit']
TS = linspace(0.25, 0.75, 11)
FASTCORE_THRESHOLDS = TS*2*log(2)
TINIT_THRESHOLDS = TS*2*log(2)

if DATASET == 'Desai-GTEx':
    THREAD_NUMBER = 5

elif DATASET == 'Blanco-Melo':
    THREAD_NUMBER = 4

PROTECTED = [['biomass_human', 'VBOF', 'EX_VBOF']] + list(UPTAKE_DRAINS['HAM Medium'])

