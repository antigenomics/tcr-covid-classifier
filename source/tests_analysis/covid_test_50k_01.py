import os
from multiprocessing import Manager, Pool

os.chdir('/home/ekvlasova/covid-t-cell-classifier')

import pandas as pd

from scipy.stats import fisher_exact
from tqdm import tqdm
from multipy.fwer import hochberg

clone_to_pval = Manager().dict()


def process_one_clone(clone):
    cmv_clone = covid_data[covid_data[clone] != 0].shape[0]
    cmv_no_clone = covid_data[covid_data[clone] == 0].shape[0]
    no_cmv_clone = healthy_data[healthy_data[clone] != 0].shape[0]
    no_cmv_no_clone = healthy_data[healthy_data[clone] == 0].shape[0]
    clone_to_pval[clone] = fisher_exact([[cmv_clone, cmv_no_clone], [no_cmv_clone, no_cmv_no_clone]])[1]
    print(f'dict size is {len(clone_to_pval)}')


if __name__ == "__main__":
    um = pd.read_csv('data/standardized_log_exp_usage_matrix_joint_new.csv').drop(columns=['Unnamed: 0']).fillna(0)

    clonotype_matrix = pd.read_csv('data/clonotype_matrix_50k.csv')
    clonotype_matrix = clonotype_matrix.drop_duplicates().set_index('cdr3aa').T.reset_index().rename(
        columns={'index': 'run'})
    clonotype_matrix = clonotype_matrix.merge(um[['run', 'covid']])
    covid_data = clonotype_matrix[clonotype_matrix.covid == 'covid']
    healthy_data = clonotype_matrix[clonotype_matrix.covid == 'healthy']
    clonotype_matrix = clonotype_matrix.drop(columns=['covid', 'run'])
    with Pool(80) as p:
        p.map(process_one_clone, clonotype_matrix.columns)

    pvals = []
    for clone in tqdm(clonotype_matrix.columns):
        pvals.append(clone_to_pval[clone])
    significant_pvals = hochberg(pvals, alpha=0.05)
    significant_clones = []
    for pval, clone in zip(significant_pvals, clonotype_matrix.columns):
        if pval:
            significant_clones.append(clone)
    pd.DataFrame(data={'clone': significant_clones}).to_csv('data/covid_significant_clones_50k.csv', index=False)
