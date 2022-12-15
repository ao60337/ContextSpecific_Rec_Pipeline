from config_variables import DATASET

MODEL_PATH = 'utils/human1_models/HumanGEM_COVID19.xml'
CONSISTENT_MODEL_PATH = 'utils/human1_models/HumanGEM_consistent_COVID19.xml'

OMICS_DATA_PATH = 'datasets/Desai-GTEx_ensembl.csv'
# OMICS_DATA_PATH = 'datasets/Blanco-Melo_ensembl.csv'

if DATASET == 'Desai-GTEx':
    TROPPO_RESULTS_PATH = 'results/troppo/Desai-GTEx'
    MODEL_RESULTS_PATH = 'results/human1_models/Desai-GTEx'

elif DATASET == 'Blanco-Melo':
    TROPPO_RESULTS_PATH = 'results/troppo/Blanco-Melo'
    MODEL_RESULTS_PATH = 'results/human1_models/Blanco-Melo'

