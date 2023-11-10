import os
from multiprocessing import Manager, Pool

os.chdir('/home/ekvlasova/covid-t-cell-classifier')

import pandas as pd

from scipy.stats import fisher_exact, binom, chi2_contingency
from tqdm import tqdm
from multipy.fwer import hochberg

clone_to_pval = Manager().dict()


def get_sum_clones_by_runs(runs=None):
    if runs is None:
        return run_to_number_of_clonotypes.number_of_clones.sum()
    return run_to_number_of_clonotypes.loc[runs, :].number_of_clones.sum()


def process_one_clone_binom(clone):
    clone, healthy_data, covid_data = clone
    covid_clone = covid_data[clone].sum()
    N = get_sum_clones_by_runs(covid_data.run)
    P = (covid_data[clone].sum() + healthy_data[clone].sum()) / get_sum_clones_by_runs()
    p_val = 1 - binom.cdf(covid_clone, N, P)
    clone_to_pval[clone] = p_val
    print(len(clone_to_pval))


def process_one_clone_fisher(clone):
    clone, healthy_data, covid_data = clone
    covid_clone = covid_data[clone].sum()
    covid_no_clone = get_sum_clones_by_runs(covid_data.run) - covid_clone
    no_covid_clone = healthy_data[clone].sum()
    no_covid_no_clone = get_sum_clones_by_runs(healthy_data.run) - no_covid_clone
    res = chi2_contingency([[covid_clone, covid_no_clone], [no_covid_clone, no_covid_no_clone]])
    if sum(map(lambda x: x < 5, res[3].flatten())) != 0:
        res = fisher_exact([[covid_clone, covid_no_clone], [no_covid_clone, no_covid_no_clone]])

    clone_to_pval[clone] = res[1]
    print(f'dict size is {len(clone_to_pval)}')


def platform_test(adaptive_data, fmba_data, clonotype_matrix, pval_save_path=None, fisher=True):
    arguments = [(x, adaptive_data, fmba_data) for x in clonotype_matrix.columns]
    print(f'Testing started')
    with Pool(40) as p:
        if fisher:
            p.map(process_one_clone_fisher, arguments)
        else:
            p.map(process_one_clone_binom, arguments)
    pvals = []
    for clone in tqdm(clonotype_matrix.columns):
        pvals.append(clone_to_pval[clone])
    significant_pvals = hochberg(pvals, alpha=0.05)
    if pval_save_path is not None:
        pd.DataFrame(data={'clone': clonotype_matrix.columns, 'pval': pvals}).to_csv(pval_save_path)
    significant_clones = []
    for pval, clone in zip(significant_pvals, clonotype_matrix.columns):
        if pval:
            significant_clones.append(clone)
    return significant_clones


def covid_test_matrix_based(clonotype_matrix, desc_path, save_path, pval_save_path=None, n=500000, platform=None,
                            fisher=True):
    desc = pd.read_csv(desc_path).drop(columns=['Unnamed: 0'])
    if platform is not None:
        desc = desc[desc.platform == platform]
    desc = desc[desc.covid == 'healthy']
    clonotype_matrix = clonotype_matrix.head(n)
    if 'cdr3aa' in clonotype_matrix.columns:
        clonotype_matrix = clonotype_matrix.drop_duplicates().set_index('cdr3aa').T.reset_index().rename(
            columns={'index': 'run'})
    if 'Unnamed: 0' in clonotype_matrix.columns:
        clonotype_matrix = clonotype_matrix.drop(columns=['Unnamed: 0'])
    clonotype_matrix = clonotype_matrix.merge(desc[['run', 'platform']])
    print(clonotype_matrix.shape)
    fmba_data = clonotype_matrix[clonotype_matrix.platform == 'fmba']
    adaptive_data = clonotype_matrix[clonotype_matrix.platform == 'adaptive']
    clonotype_matrix = clonotype_matrix.drop(columns=['platform', 'run'])
    print(f'Finished preparing data. Adaptive data shape: {adaptive_data.shape}, FMBA: {fmba_data.shape}')

    significant_clones = platform_test(adaptive_data, fmba_data, clonotype_matrix, pval_save_path, fisher)

    pd.DataFrame(data={'clone': significant_clones}).to_csv(save_path, index=False)


def platform_test_for_beta_chain_adaptive(n=500000):
    clm = pd.read_csv('data/adaptive_clone_results/clonotype_matrix_fmba_adaptive_top_500k_1_mismatch.csv')
    print(f'Finished reading clonotype matrix')
    covid_test_matrix_based(
        clonotype_matrix=clm,
        desc_path=f'data/standardized_log_exp_usage_matrix_joint_new.csv',
        save_path=f'data/anomaly_clones/platform_clones_healthy_top_500k_1_mismatch_fisher.csv',
        pval_save_path=f'data/anomaly_clones/platform_pvals_healthy_top_500k_1_mismatch_fisher.csv',
        n=n,
        fisher=True,
    )


if __name__ == "__main__":
    run_to_number_of_clonotypes = pd.read_csv('data/run_to_number_of_clonotypes.csv').set_index('run')
    platform_test_for_beta_chain_adaptive()
