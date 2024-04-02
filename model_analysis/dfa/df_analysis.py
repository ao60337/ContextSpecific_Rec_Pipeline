import traceback

import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api
from cobra.sampling import OptGPSampler
from scipy.stats import hypergeom, ks_2samp


def bootstrapCI(rxn):
    """
    Calculate the confidence interval of a reaction

    Parameters
    ----------
    rxn

    Returns
    -------
    int: 1 if the reaction is significantly different from the mean, 0 otherwise.

    """
    bsci = []

    for i in range(1000):
        bt_samp = rxn.sample(1000, replace=True)
        bsci.append(bt_samp.mean())

    ci_low = np.percentile(bsci, 2.5)
    ci_high = np.percentile(bsci, 97.5)

    if ci_low > 0 or ci_high < 0:
        return 1
    else:
        return 0


def kstest(samples_healthy: pd.DataFrame, samples_infected: pd.DataFrame, dataset_name: str, target_tissue: str):
    """
    Calculate the K-S test to detect significantly altered reactions fluxes.
    Results are saved in a csv file.

    Parameters
    ----------
    samples_healthy: pd.DataFrame
        The samples of the healthy tissue.
    samples_infected: pd.DataFrame
        The samples of the infected tissue.
    dataset_name: str
        The name of the dataset.
    target_tissue: str
        The name of the target tissue.

    Returns
    -------
    pd.DataFrame: The results of the K-S test for each reaction.

    """
    rxns1 = set(samples_infected.columns)
    rxns2 = set(samples_healthy.columns)

    rxn_c = rxns1.intersection(rxns2)

    pvals = []
    rxnid = []
    fc = []

    for rxn in rxn_c:
        data1 = samples_infected[rxn].round(decimals=4)
        data2 = samples_healthy[rxn].round(decimals=4)

        data1 = data1.sample(n=1000)
        data2 = data2.sample(n=1000)

        if (data1.std() != 0 and data1.mean() != 0) or (data2.std() != 0 and data2.mean() != 0):
            kstat, pval = ks_2samp(data1, data2)

            foldc = (data1.mean() - data2.mean()) / abs(data1.mean() + data2.mean())

            pvals.append(pval)
            rxnid.append(rxn)
            fc.append(foldc)

    data_mwu = pd.DataFrame({'Reaction': rxnid, 'Pvalue': pvals})
    data_mwu = data_mwu.set_index('Reaction')

    reject, padj, _, _ = statsmodels.stats.multitest.multipletests(data_mwu['Pvalue'], alpha=0.05, method='fdr_bh',
                                                                   is_sorted=False, returnsorted=False)

    data_mwu['Padj'] = padj
    data_mwu['Reject'] = reject
    data_mwu['FC'] = fc

    data_sigFC = data_mwu.loc[(abs(data_mwu['FC']) > 0.82) & (data_mwu['Padj'] < 0.05), :]

    rxns1 = set(samples_infected.columns)
    rxns2 = set(samples_healthy.columns)

    rxn_in1 = rxns1.difference(rxns2)
    rxn_in2 = rxns2.difference(rxns1)

    act = []
    rep = []

    for rx in rxn_in1:  # Activated reactions
        sig = bootstrapCI(samples_infected[rx])

        if sig == 1:
            act.append(rx)

    for rx in rxn_in2:  # Activated reactions
        sig = bootstrapCI(samples_healthy[rx])

        if sig == 1:
            rep.append(rx)

    df_abs = pd.DataFrame({'Reaction': act + rep, 'Padj': np.zeros(len(act + rep))})
    df_abs = df_abs.set_index('Reaction')
    data_return = data_sigFC + df_abs

    file = 'Differential_Flux_Analysis/%s/%s/%s_DFA_reaction_result2.csv' \
           % (dataset_name, target_tissue, target_tissue)
    data_return.to_csv(file)

    return data_return.index.to_list()


def pathway_enrichment(rxnlist: list, dataset_name: str, target_tissue: str):
    """
    Maps significantly altered reactions to pathways using the subsystems from the HumanGEM model.
    Results are saved in csv and jpg files.

    Parameters
    ----------
    rxnlist: list
        The list of reactions to be mapped to pathways.
    dataset_name: str
        The name of the dataset.
    target_tissue: str
        The name of the target tissue.

    """
    subs = pd.read_csv('HumanGEM_Subsystems_reduced.csv', sep=',')
    dataset = pd.DataFrame()

    for path in subs['Subsystem'].unique():
        reaction_set = subs.loc[subs['Subsystem'] == path, 'Reaction']

        rxn = reaction_set.reset_index(drop=True)

        df_temp = pd.DataFrame({path: rxn})

        dataset = pd.concat([dataset, df_temp], axis=1)

    listrxnSize = []
    setSize = []

    d = [g for g in rxnlist]

    for col in dataset.columns:
        df = pd.DataFrame({'Reaction': dataset[col]})

        out = []

        for reac in df['Reaction']:
            if reac in rxnlist:
                out.append(reac)
                d.remove(reac)

        listrxnSize.append(len(out))
        setSize.append(len(dataset[col].dropna()))

    hyperdata = pd.DataFrame({'Pathways': dataset.columns, 'ListReactions': listrxnSize, 'SetSize': setSize})

    hits = hyperdata['ListReactions']
    pool = hyperdata['SetSize']

    allrxns = hyperdata['SetSize'].sum()
    targetrxns = hyperdata['ListReactions'].sum()

    pvalList = []

    for h, p in zip(hits, pool):
        rv = hypergeom(allrxns - p, p, targetrxns)

        pval = rv.pmf(h)

        pvalList.append(pval)

    hyperdata['P-value'] = pvalList

    reject, padj, _, _ = statsmodels.stats.multitest.multipletests(hyperdata['P-value'], alpha=0.05, method='fdr_bh',
                                                                   is_sorted=False, returnsorted=False)

    hyperdata['P-value_adj'] = padj
    hyperdata['Reject'] = reject

    hyperdata_sig = hyperdata[(hyperdata['Reject']) & (hyperdata['ListReactions'] != 0)]

    # hyperdata_sorted = hyperdata_sig.sort_values(by='P-value_adj', ascending=False)
    hyperdata_sorted = hyperdata_sig.sort_values(by='ListReactions', ascending=True)
    hyperdata_sorted = hyperdata_sorted[hyperdata_sorted['Pathways'] != 'Transport reactions']
    hyperdata_sorted = hyperdata_sorted[hyperdata_sorted['Pathways'] != 'Exchange/demand reactions']
    hyperdata_sorted = hyperdata_sorted[hyperdata_sorted['Pathways'] != 'Pool reactions']
    hyperdata_sorted = hyperdata_sorted[hyperdata_sorted['Pathways'] != 'Isolated']
    hyperdata_sorted = hyperdata_sorted[hyperdata_sorted['Pathways'] != 'Drug metabolism']

    plt.figure(figsize=(15, 15))

    sc = plt.scatter(hyperdata_sorted['P-value_adj'], np.arange(0, len(hyperdata_sorted['Pathways'])),
                     s=hyperdata_sorted['ListReactions'], color=(0.9, 0.3, 0.1, 0.9))

    plt.xlabel('Adjusted p-value')

    plt.yticks(np.arange(0, len(hyperdata_sorted['Pathways'])), labels=hyperdata_sorted['Pathways'])

    handles, labels = sc.legend_elements(prop="sizes", alpha=0.6)

    plt.legend(handles, labels, bbox_to_anchor=(1.6, 1.02), loc='upper right', title="Reactions")

    plt.tight_layout()

    plt.savefig('results/figures/%s_DFA_result2.jpg' % target_tissue, dpi=600)

    hyperdata_sorted.to_csv('results/csv/%s_DFA_pathway_result.csv' % target_tissue)


def dfa(data: str, target_tissue: str, tissue_models: dict):
    """
    Performs differential flux analysis to compare the flux  of healthy and infected models.

    Parameters
    ----------
    data: str
        The name of the dataset.
    target_tissue: str
        The name of the target tissue.
    tissue_models: dict
        The dictionary with the model names for that target tissue.

    """
    model_healthy = cobra.io.read_sbml_model(f'../models/{tissue_models["healthy"]}.csv')
    model_infected = cobra.io.read_sbml_model(f'../models/{tissue_models["infected"]}.csv')
    #
    model_healthy.objective = 'biomass_human'
    model_infected.objective = 'VBOF'
    #
    print('Starting ACHRSampler Healthy')
    achr_healthy = OptGPSampler(model_healthy, processes=12, thinning=100)

    print('Starting ACHRSampler Infected')
    achr_infected = OptGPSampler(model_infected, processes=12, thinning=100)
    #
    print('Starting Sampling Healthy')
    samples_healthy = achr_healthy.sample(1000)
    print('Starting Sampling Infected')
    samples_infected = achr_infected.sample(1000)

    samples_healthy.to_csv('results/csv/%s_samples.csv'
                           % (tissue_models['healthy']))
    samples_infected.to_csv('results/csv/%s_samples.csv'
                            % (tissue_models['infected']))

    print('Starting KS-test')
    ks_result = kstest(samples_healthy, samples_infected, data, target_tissue)

    print('Starting Pathway Enrichment')
    pathway_enrichment(list(ks_result['Reaction']), data, target_tissue)

    print('Finished')


if __name__ == '__main__':
    tissues_desai = {
        'Lung': {'healthy': 'HumanGEM_Lung_Healthy_local2_1_3_2_fastcore_t0',
                 'infected': 'HumanGEM_Lung_COVID19_local2_1_3_2_fastcore_t0'},
        'Heart': {'healthy': 'HumanGEM_Heart_Healthy_local2_1_3_2_fastcore_t0',
                  'infected': 'HumanGEM_Heart_COVID19_local2_1_3_2_fastcore_t0'},
        'Intestine': {'healthy': 'HumanGEM_Intestine_Healthy_local2_1_3_2_fastcore_t0',
                      'infected': 'HumanGEM_Intestine_COVID19_local2_1_3_2_fastcore_t0'},
        'Kidney': {'healthy': 'HumanGEM_Kidney_Healthy_local2_1_3_2_fastcore_t0',
                   'infected': 'HumanGEM_Kidney_COVID19_local2_1_3_2_fastcore_t0'},
        'Liver': {'healthy': 'HumanGEM_Liver_Healthy_local2_1_3_2_fastcore_t0',
                  'infected': 'HumanGEM_Liver_COVID19_local2_1_3_2_fastcore_t0'}
    }

    for tissue in tissues_desai.keys():
        print('-----------------------------------------')
        print('DFA for %s models' % tissue)
        print('-----------------------------------------')

        # noinspection PyBroadException
        try:
            dfa(data='Desai-GTEx', target_tissue=tissue, tissue_models=tissues_desai[tissue])
        except:
            print('DFA for %s failed' % tissue)
            traceback.print_exc()
            # continue