import traceback
import warnings
import cobra.flux_analysis
import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions
from cobamp.wrappers import KShortestMCSEnumeratorWrapper, KShortestGeneticMCSEnumeratorWrapper

warnings.simplefilter('ignore')


def read_tissue_models(tissue_models: dict) -> object:
    """
    Opens the models for a given tissue.

    Parameters
    ----------
    tissue_models: dict
        The dictionary with the model names for that target tissue.

    Returns
    -------
    model_infected: cobra.Model
        The model for the infected tissue.

    """
    model_infected = read_sbml_model('../models/%s.xml' % (tissue_models['infected']))
    model_infected.remove_reactions(find_blocked_reactions(model_infected))
    model_infected.objective = 'VBOF'

    print('Model loaded and blocked reactions removed.')

    return model_infected


def get_wt_simulation(model_infected: cobra.Model):
    """
    Calculates the wild-type simulation for the healthy and infected tissues.

    Parameters
    ----------
    model_infected: cobra.Model
        The model for the infected tissue.

    Returns
    -------
    wt_infected: float
        The wild-type simulation for the infected tissue.

    """
    wt_infected = model_infected.optimize().objective_value

    print('Wild-type simulations for healthy and infected models calculated.')

    return wt_infected


def get_non_targets(cobra_model: cobra.Model, wt_simulation: float,
                    model_objective: str, model_name: str):
    """
    Gets the non-targets for the model.

    Parameters
    ----------
    cobra_model: cobra.Model
        The model for the tissue.
    wt_simulation: float
        The wild-type simulation for the tissue.
    model_objective: str
        The objective reaction of the model.
    model_name: str
        The name of the model.

    Returns
    -------
    target_flux_space: dict
        The target flux space for the model.
    non_targets: list
        The non-targets for the model.

    """
    target_flux_space = {k: (-v, None) for k, v in cobra_model.medium.items()}
    target_flux_space.update({model_objective: (0.001 * wt_simulation, None)})

    non_targets = list(target_flux_space.keys())

    biomass = ['HMR_10062', 'HMR_10063', 'HMR_10064', 'HMR_10065', 'biomass_human', 'VBOF', 'EX_VBOF']

    reaction_deletions = cobra.flux_analysis.single_reaction_deletion(cobra_model, method='pfba', processes=4)

    essential_reactions = []
    for i in range(reaction_deletions.shape[0]):
        if round(reaction_deletions.iloc[i]['growth'], 6) == 0.0:
            essential_reactions.append(str(reaction_deletions.iloc[i]['ids']).replace("{'", '').replace("'}", ''))

    essentials = pd.DataFrame.from_dict({'Essential_Reactions': essential_reactions})
    essentials.to_csv('results/%s_essential_reactions.csv' % model_name)

    for reaction in cobra_model.reactions:
        if len(reaction.compartments) == 2 or reaction.id in biomass or reaction.id in essential_reactions:
            non_targets.append(reaction.id)

    for drain in cobra_model.exchanges:
        if drain.id not in non_targets:
            non_targets.append(drain.id)

    print('Target flux space and non-target reactions calculated.')

    return target_flux_space, non_targets


def find_minimal_cut_sets(model: cobra.Model, target_flux_space,
                          non_targets, model_name: str):
    """
    Finds the minimal cut sets for the model using Cobamp.

    Parameters
    ----------
    model: cobra.Model
        The model for the tissue.
    target_flux_space: dict
        The target flux space for the model.
    non_targets: list
        The non-targets for the model.
    model_name: str
        The name of the model.

    Returns
    -------
    cut_set_list: list
        The minimal cut sets for the model.

    """
    iterations = 2  # how many sizes / solutions to generate
    by_size = True  # generate solutions one size or one set of KOs at a time?

    mcs = KShortestMCSEnumeratorWrapper(model=model,
                                        algorithm_type='kse_populate' if by_size else 'kse_iterative',
                                        stop_criteria=iterations,
                                        big_m=False,
                                        target_yield_space_dict={},
                                        target_flux_space_dict=target_flux_space,
                                        excluded_solutions=[[k] for k in non_targets],
                                        solver='CPLEX')

    enumerator = mcs.get_enumerator()

    solutions = [result for result in enumerator]
    mcs = [set(solution.keys()) for solution in solutions[1]]

    cut_set_list = []

    for cut_set in mcs:
        cut_set_list.append(list(cut_set))

    cs_df = pd.DataFrame(cut_set_list, columns=['Reaction 1', 'Reaction2'])
    cs_df.to_csv('results/%s_msc.csv' % model_name)

    print('Minimal cut sets for the %s model calculated.' % model_name)

    return cut_set_list


def mcs_pipeline(dataset_name: str, tissue_name: str, models: dict):
    """
    The pipeline for the MCS analysis.
    Outputs a csv file that contains the minimal cut sets for the infected tissue not in the healthy tissues.

    Parameters
    ----------
    dataset_name: str
        The name of the dataset.
    tissue_name: str
        The name of the tissue.
    models: dict
        The models for the tissue.

    """
    result_path = 'Drug_Targeting/%s/%s' % (dataset_name, tissue_name)

    infected_tissue = read_tissue_models(models)
    wt_infected = get_wt_simulation(infected_tissue)

    infected_flux_space, infected_non_targets = get_non_targets(infected_tissue, wt_infected, 'VBOF',
                                                                result_path)

    infected_mcs = find_minimal_cut_sets(infected_tissue, infected_flux_space, infected_non_targets, models['infected'])

    return infected_mcs


if __name__ == '__main__':
    tissues_desai = {'Lung': {'infected': 'HumanGEM_Lung_COVID19_local2_1_3_3_fastcore_t0'},
                     'Heart': {'infected': 'HumanGEM_Heart_COVID19_local2_1_3_3_fastcore_t0'},
                     'Intestine': {'infected': 'HumanGEM_Intestine_COVID19_local2_1_3_3_fastcore_t0'},
                     'Kidney': {'infected': 'HumanGEM_Kidney_COVID19_local2_1_3_3_fastcore_t0'},
                     'Liver': {'infected': 'HumanGEM_Liver_COVID19_local2_1_3_3_fastcore_t0'}
                     }

    for tissue in tissues_desai.keys():
        print('-----------------------------------------')
        print('MCS analysis for %s models' % tissue)
        print('-----------------------------------------')

        # noinspection PyBroadException
        try:
            mcs_pipeline(dataset_name='Desai-GTEx', tissue_name=tissue, models=tissues_desai[tissue])
        except:
            print('Minimal Cut Set Analysis for %s failed' % tissue)
            traceback.print_exc()
            continue
