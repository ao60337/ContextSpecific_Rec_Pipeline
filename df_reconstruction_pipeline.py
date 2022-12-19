import sys
import os
import warnings
import traceback

# Warnings and prints are blocked to avoid massive outputs that appear after running COBAMP.
warnings.simplefilter("ignore")


# Disable
def block_print():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enable_print():
    sys.stdout = sys.__stdout__


block_print()

import cobra
import pandas as pd
import re

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import write_sbml_model
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from cobamp.utilities.parallel import batch_run
from utils.pipeline_paths import *
from utils.medium_variables import *
from utils.config_variables import *


def print_model_details(cobra_model):
    """
    Function to print the details of the currently loaded COBRA model.

    Parameters
    ----------
    cobra_model : cobra.Model

    """
    enable_print()

    transporters = []

    for reac in cobra_model.reactions:
        if len(reac.compartments) == 2:
            transporters.append(reac.id)

    print('Total Reactions:', len(cobra_model.reactions))
    print('Reactions:', (len(cobra_model.reactions)) - len(transporters) - len(cobra_model.exchanges))
    print('Transporters:', len(transporters))
    print('Exchanges:', len(cobra_model.exchanges))

    block_print()


def load_model(model_path: str, consistent_model_path: str):
    """
    This function is used to load the model.

    Parameters
    ----------
    model_path : str
        The path to the model.
    consistent_model_path : str
        The path to the model without blocked reactions.

    Returns
    -------
    model : cobra.Model
        The loaded model.
    """

    if os.path.exists(consistent_model_path):
        model = cobra.io.read_sbml_model(consistent_model_path)

    else:
        model = cobra.io.read_sbml_model(model_path)
        model.remove_reactions(find_blocked_reactions(model))
        write_sbml_model(model, consistent_model_path)

    print_model_details(model)

    for reaction_id, bound in MEDIUM_CONDITIONS.items():
        if reaction_id in model.reactions:
            model.reactions.get_by_id(reaction_id).bounds = bound

    return model


def reconstruction_function(omics_container, parameters: dict):
    """
    This function is used to run the reconstruction algorithm.

    Parameters
    ----------
    omics_container : pandas.DataFrame
        The omics data set.
    parameters : dict
        The parameters to be used for the reconstruction algorithm.

    Returns
    ----------
    rec_wrapper : Reconstruction Wrapper object with model and omics data.
    """

    def integration_fx(data_map):
        return [[k for k, v in data_map.get_scores().items() if (v is not None and v > threshold) or k in PROTECTED]]

    def score_apply(data_map):
        dm = {k: 0 if v is None else (min(v, 10) - threshold) if k not in PROTECTED else 20
              for k, v in data_map.items()}
        return dm

    threshold, rec_wrapper, method = [parameters[parameter] for parameter in
                                      ['threshold', 'reconstruction_wrapper', 'algorithm']]

    # noinspection PyBroadException
    try:
        if method == 'fastcore':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method, and_or_funcs=(min, sum),
                                              integration_strategy=('custom', [integration_fx]), solver='CPLEX')

        elif method == 'tinit':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method, and_or_funcs=(min, sum),
                                              integration_strategy=('continuous', score_apply), solver='CPLEX')

    except:
        traceback.print_exc()

        return {r: False for r in rec_wrapper.model_reader.r_ids}


def troppo_integration(model: cobra.Model, algorithm: str, threshold: float, thread_number: int):
    """
    This function is used to run the Troppo's integration algorithms.

    Parameters
    ----------
    model: cobra.Model
        The COBRA model.
    algorithm: str
        The algorithm to be used.
    threshold: float
        The threshold to be used.
    thread_number: int
        The number of threads to be used.

    Returns
    -------
    integration_results: dict
        Dataframe containing the results of the omics integration.
        Each sample as a dictionary containing a boolean value for each reaction.

    """

    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    details = [DATASET, algorithm, threshold]

    template = model.copy()

    omics_dataset = pd.read_csv(OMICS_DATA_PATH, index_col=0)

    troppo_results_path = os.path.join(TROPPO_RESULTS_PATH, 'HumanGEM_%s_%s_%s.csv') % (DATASET, algorithm,
                                                                                        round(threshold, 2))

    omics_data = TabularReader(path_or_df=omics_dataset, nomenclature='ensemble_gene_id',
                               omics_type='transcriptomics').to_containers()

    enable_print()
    print('Tabular Reader Finished.')
    block_print()

    reconstruction_wrapper = ReconstructionWrapper(template, ttg_ratio=9999,
                                                   gpr_gene_parse_function=replace_alt_transcripts)

    enable_print()
    print('Reconstruction Wrapper Finished.')
    block_print()

    parameters = {'threshold': threshold, 'reconstruction_wrapper': reconstruction_wrapper, 'algorithm': algorithm}

    batch_fastcore_res = batch_run(reconstruction_function, omics_data, parameters, threads=thread_number)

    result_dict = dict(zip([sample.condition for sample in omics_data], batch_fastcore_res))

    integration_result = pd.DataFrame.from_dict(result_dict, orient='index')
    integration_result.to_csv(troppo_results_path)

    enable_print()
    print('Omics Integration with %s (Threshold = %s) Finished.' % (details[1], details[2]))
    block_print()

    return result_dict


def sbml_model_reconstruction(model_template: cobra.Model, sample: str, integration_result_dict: dict):
    """
    This function is used to reconstruct the model based on the integration results.

    Parameters
    ----------
    model_template: cobra.Model
        The COBRA model template.
    sample: str
        The sample name.
    integration_result_dict: dict
        The integration results.
    """

    with model_template as temp_model:
        if 'COVID19' in sample:
            temp_model.objective = 'VBOF'
        else:
            temp_model.objective = 'biomass_human'

        reactions_to_deactivate = [reaction for reaction, value in
                                   integration_result_dict[sample].items() if value is False]

        for reaction in reactions_to_deactivate:
            temp_model.remove_reactions([reaction])

        enable_print()

        for reaction_id, bound in MEDIUM_CONDITIONS['HAM Medium'].items():
            if reaction_id in model_template.reactions:
                model_template.reactions.get_by_id(reaction_id).bounds = bound
            else:
                print(reaction_id, 'exchange not found in the model.')

        model_name = sample.split('_')[3] + '/' + sample.split('_')[1] + '/' + sample + '.xml'

        cobra.io.write_sbml_model(temp_model, os.path.join(MODEL_RESULTS_PATH, model_name))

        print('Model reconstruction for %s finished.' % sample)
        print_model_details(temp_model)

        block_print()


def default_reconstruction_pipeline():
    """
    This function is used to run the reconstruction pipeline.

    In the end this function generates a SBML file of the tissue-specific models that resulted from the omics data
    integration with Troppo.

    All the parameters required to run this pipeline can be defined in the pipeline_paths.py, config_variables.py,
    and medium_variables.py files.
    """
    enable_print()

    print('-------------------------------------------------------------------------------------------------------')
    print('--------------------------------------- Loading template model. ---------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    block_print()

    template_model = load_model(MODEL_PATH, CONSISTENT_MODEL_PATH)

    enable_print()

    print('-------------------------------------------------------------------------------------------------------')
    print('------------------------------- Starting Omics Integration with Troppo. -------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    block_print()

    integration_result = {}

    for algorithm in ALGORITHMS:

        enable_print()

        print('Omics Integration with %s Started.' % algorithm)
        print('-------------------------------------------------------------------------------------------------------')

        block_print()

        if algorithm == 'fastcore':
            for threshold in FASTCORE_THRESHOLDS:
                troppo_result = troppo_integration(template_model, algorithm, threshold, THREAD_NUMBER_FASTCORE)
                for sample in list(troppo_result.keys()):
                    integration_result['HumanGEM_%s_%s_%s' % (sample, algorithm, round(threshold, 2))] \
                        = troppo_result[sample]

            enable_print()

            print('----------------------------------------------------------'
                  '---------------------------------------------')

            block_print()

        elif algorithm == 'tinit':
            for threshold in TINIT_THRESHOLDS:
                troppo_result = troppo_integration(template_model, algorithm, threshold, THREAD_NUMBER_TINIT)
                for sample in list(troppo_result.keys()):
                    integration_result['HumanGEM_%s_%s_%s' % (sample, algorithm, round(threshold, 2))] = \
                        troppo_result[sample]

            enable_print()

            print('----------------------------------------------------------'
                  '---------------------------------------------')

            block_print()

    enable_print()

    print('-------------------------------------------------------------------------------------------------------')
    print('----------------------- Starting Reconstruction of the Context-Specific Models. -----------------------')
    print('-------------------------------------------------------------------------------------------------------')

    block_print()

    for sample_name in list(integration_result.keys()):
        enable_print()

        print('Context-specific model reconstruction for %s started.' % sample_name)
        print('-------------------------------------------------------------------------------------------------------')

        sbml_model_reconstruction(template_model, sample_name, integration_result)

        print('-------------------------------------------------------------------------------------------------------')

        block_print()

    # TODO: Add the implementation of the local threshold functions to the pipeline.
    # TODO: Add the option for more integration algorithms in the pipeline.
    # TODO: Add a function to implement gap-filling to the human1_models generated with troppo.
    # TODO: Add a function to evaluate task performance of the human1_models generated with troppo.

    enable_print()

    print('-------------------------------------------------------------------------------------------------------')
    print('------------------------------------------ Pipeline Finished ------------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    block_print()

    return integration_result


if __name__ == '__main__':
    default_reconstruction_pipeline()
