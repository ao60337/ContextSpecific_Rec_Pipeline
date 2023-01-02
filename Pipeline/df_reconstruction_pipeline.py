import os
import pandas as pd

from utils.pipeline_paths import *
from utils.config_variables import *
from model_handle import load_model, sbml_model_reconstruction
from troppo_integration import troppo_integration


def default_reconstruction_pipeline():
    """
    This function is used to run the reconstruction pipeline.

    In the end this function generates a SBML file of the tissue-specific reconstructed_models that resulted from the
    omics data integration with Troppo.

    All the parameters required to run this pipeline can be defined in the ***pipeline_paths.py***,
    ***config_variables.py***, and ***medium_variables.py files***.
    """
    print('-------------------------------------------------------------------------------------------------------')
    print('--------------------------------------- Loading template model. ---------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    template_model = load_model(MODEL_PATH, CONSISTENT_MODEL_PATH)

    print('-------------------------------------------------------------------------------------------------------')
    print('------------------------------- Starting Omics Integration with Troppo. -------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    omics_data = pd.read_csv(OMICS_DATA_PATH, index_col=0)

    integration_result = {}

    for algorithm in ALGORITHMS:

        print('Omics Integration with %s Started.' % algorithm)
        print('-------------------------------------------------------------------------------------------------------')

        if algorithm == 'fastcore':
            thred_number = THREAD_NUMBER_FASTCORE

        elif algorithm == 'tinit':
            thred_number = THREAD_NUMBER_TINIT

        else:
            return 'Algorithm not supported by this pipeline.'

        for threshold in INTEGRATION_THRESHOLDS:
            troppo_result = troppo_integration(template_model, algorithm, threshold, thred_number, omics_data)

            for sample in list(troppo_result.keys()):
                th = str(round(threshold, 2)).replace('.', '_')
                integration_result['%s_%s_%s_%s_t%s' % (MODEL, sample, algorithm, THRESHOLDING_STRATEGY,
                                                        th)] = troppo_result[sample]

            print('----------------------------------------------------------'
                  '---------------------------------------------')

    troppo_results_path = os.path.join(TROPPO_RESULTS_PATH, '%s_%s_%s.csv') % (MODEL, DATASET, THRESHOLDING_STRATEGY)
    integration_dataframe = pd.DataFrame.from_dict(integration_result, orient='index')
    integration_dataframe.to_csv(troppo_results_path)

    print('----------------------- Starting Reconstruction of the Context-Specific Models. -----------------------')
    print('-------------------------------------------------------------------------------------------------------')

    for sample_name in list(integration_result.keys()):
        print('Context-specific model reconstruction for %s started.' % sample_name)
        print('-------------------------------------------------------------------------------------------------------')

        sbml_model_reconstruction(template_model, sample_name, integration_result)

        print('-------------------------------------------------------------------------------------------------------')

    # TODO: Add the implementation of the local threshold functions to the pipeline.
    # TODO: Add the option for more integration algorithms in the pipeline.
    # TODO: Add a function to implement gap-filling to the reconstructed_models generated with troppo.
    # TODO: Add a function to evaluate task performance of the reconstructed_models generated with troppo.

    print('------------------------------------------ Pipeline Finished ------------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    return integration_result


if __name__ == '__main__':
    default_reconstruction_pipeline()
