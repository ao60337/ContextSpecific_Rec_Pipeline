import cobra
import os

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import write_sbml_model
from utils.medium_variables import MEDIUM_CONDITIONS
from utils.config_variables import OBJECTIVE, MEDIUM_NAME
from utils.pipeline_paths import MODEL_RESULTS_PATH


def print_model_details(cobra_model):
    """
    Function to print the details of the currently loaded COBRA model.

    Parameters
    ----------
    cobra_model : cobra.Model

    """
    transporters = []

    for reac in cobra_model.reactions:
        if len(reac.compartments) == 2:
            transporters.append(reac.id)

    print('Total Reactions:', len(cobra_model.reactions))
    print('Reactions:', (len(cobra_model.reactions)) - len(transporters) - len(cobra_model.exchanges))
    print('Transporters:', len(transporters))
    print('Exchanges:', len(cobra_model.exchanges))


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
        temp_model.objective = OBJECTIVE

        reactions_to_deactivate = [reaction for reaction, value in
                                   integration_result_dict[sample].items() if value is False]

        for reaction in reactions_to_deactivate:
            temp_model.remove_reactions([reaction])

        for reaction_id, bound in MEDIUM_CONDITIONS[MEDIUM_NAME].items():
            if reaction_id in temp_model.reactions:
                temp_model.reactions.get_by_id(reaction_id).bounds = bound
            else:
                print(reaction_id, 'exchange not found in the model.')

        model_name = sample.split('_')[1] + '/' + sample + '.xml'

        cobra.io.write_sbml_model(temp_model, os.path.join(MODEL_RESULTS_PATH, model_name))

        print('Model reconstruction for %s finished.' % sample)
        print_model_details(temp_model)

