import os
from multiprocessing import Manager, Pool

os.chdir('/home/ekvlasova/covid-t-cell-classifier')

import pandas as pd

from scipy.stats import fisher_exact
from tqdm import tqdm
from multipy.fwer import hochberg

clone_to_pval = Manager().dict()


def get_sum_clones_by_runs(runs):
    return run_to_number_of_clonotypes.loc[runs, :].number_of_clones.sum()


def process_one_clone(clone):
    cmv_clone = cmv_data[clone].sum()
    cmv_no_clone = get_sum_clones_by_runs(cmv_data.run) - cmv_clone
    no_cmv_clone = healthy_data[clone].sum()
    no_cmv_no_clone = get_sum_clones_by_runs(healthy_data.run) - cmv_clone
    clone_to_pval[clone] = fisher_exact([[cmv_clone, cmv_no_clone], [no_cmv_clone, no_cmv_no_clone]])[1]
    print(f'dict size is {len(clone_to_pval)}')


if __name__ == "__main__":
    desc_hip = pd.read_csv('/projects/fmba_covid/1_data_links/hip_full_prep.txt', sep='\t').rename(
        columns={'file_name': 'run'})
    run_to_number_of_clonotypes = pd.read_csv('data/run_to_number_of_clonotypes.csv')
    run_to_number_of_clonotypes = run_to_number_of_clonotypes[(run_to_number_of_clonotypes.run.str.contains('HIP')) | (
        run_to_number_of_clonotypes.run.str.contains('Keck'))].set_index('run')
    clonotype_matrix = pd.read_csv('data/clonotype_matrix_hip_10k_top_1_mismatch_uniform.csv')
    # clonotype_matrix = pd.read_csv('data/clonotype_matrix_50k_1_mismatch_roulette.csv')
    clonotype_matrix = clonotype_matrix[['cdr3aa'] + [x for x in clonotype_matrix.columns if 'Keck' in x or 'HIP' in x]]
    clonotype_matrix = clonotype_matrix.drop_duplicates().set_index('cdr3aa').T.reset_index().rename(
        columns={'index': 'run'})
    clonotype_matrix = clonotype_matrix.merge(desc_hip[['run', 'cmv']])

    cmv_data = clonotype_matrix[clonotype_matrix.cmv == '+']
    healthy_data = clonotype_matrix[clonotype_matrix.cmv == '-']
    clonotype_matrix = clonotype_matrix.drop(columns=['cmv', 'run'])
    with Pool(48) as p:
        p.map(process_one_clone, clonotype_matrix.columns)

    pvals = []
    for clone in tqdm(clonotype_matrix.columns):
        pvals.append(clone_to_pval[clone])
    pd.DataFrame(data={'clone': clonotype_matrix.columns, 'pvals': pvals}).to_csv(
        'data/cmv_clones_pval_uniform_10k_1mismatch.csv', index=False)
    significant_pvals = hochberg(pvals, alpha=0.05)
    significant_clones = []
    for pval, clone in zip(significant_pvals, clonotype_matrix.columns):
        if pval:
            significant_clones.append(clone)
    pd.DataFrame(data={'clone': significant_clones}).to_csv(
        'data/cmv_significant_clones_10k_hip_uniform_1_mismatch.csv', index=False)
