from utils.config_variables import PROTECTED
from utils.pipeline_paths import MODEL_TASKS_PATH, MODEL_TASKS_RESULTS_PATH
from json import JSONEncoder
from troppo.tasks.core import TaskEvaluator
from troppo.tasks.task_io import JSONTaskIO
from cobamp.utilities.parallel import batch_run
import cobra
import sys
import os
import warnings

warnings.simplefilter("ignore")


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


blockPrint()


def read_task_file(generic_model: cobra.Model) -> list:
    """
    This function reads the task file and returns a list of tasks.

    Parameters
    ----------
    generic_model: cobra.Model
        The generic model to be used for the task evaluation.

    Returns
    -------
    list: list of tasks

    """
    return [task for task in JSONTaskIO().read_task(MODEL_TASKS_PATH)
            if len((set(task.inflow_dict) | set(task.outflow_dict)) -
                   set([metabolite.id for metabolite in generic_model.metabolites])) == 0]


def write_result_jason(model_name: str, task_evaluation_result: dict):
    """
    Function to write the task evaluation results to a json file.

    Parameters
    ----------
    model_name: str
        The name of the model that was evaluated.
    task_evaluation_result: dict
        The result of the task evaluation.

    """
    model_data = model_name.split('_')
    model_path = f'{model_data[1]}/{model_name}.json'
    result_path = os.path.join(MODEL_TASKS_RESULTS_PATH, model_path)

    with open(result_path, 'w') as f:
        f.write(JSONEncoder().encode([(k, v) for k, v in task_evaluation_result.items()]))


def task_eval(model_template: cobra.Model, integration_results: dict):
    """
    This function performs task evaluation on the context-specific models that were generated using troppo.

    Parameters
    ----------
    model_template: cobra.Model
        The generic model to be used for the task evaluation.
    integration_results: dict
        The results of the integration process.

    """
    task_list = read_task_file(model_template)

    for task in task_list:
        task.inflow_dict = {k: v if k not in task.outflow_dict.keys() else [-1000, 1000] for k, v in
                            task.inflow_dict.items()}

        task.outflow_dict = {k: v for k, v in task.outflow_dict.items() if k not in task.inflow_dict.items()}

    for task in task_list:
        task.mandatory_activity = []

    for k in model_template.boundary:
        k.knock_out()

    all_reactions = set([reaction.id for reaction in model_template.reactions])

    for model_name, result in integration_results.items():
        with model_template as context_specific_model:
            active_reactions = set([k for k, v in result.items() if v])
            protected_reactions = active_reactions + PROTECTED
            reactions_to_remove = all_reactions - protected_reactions

            for reaction_to_remove in reactions_to_remove:
                context_specific_model.reactions.get_by_id(reaction_to_remove).knock_out()

            evaluation = TaskEvaluator(model=context_specific_model, tasks=task_list, solver='CPLEX')
            task_names = evaluation.tasks
            batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': evaluation}, threads=1)

        task_result = {model_name: (v[0], v[2]) for model_name, v in dict(zip(task_names, batch_res_tasks)).items()}

        enablePrint()
        print(model_name, len(protected_reactions),
              len([v for model_name, v in task_result.items() if v[0]]),
              'tasks completed.')

        write_result_jason(model_name, task_result)




