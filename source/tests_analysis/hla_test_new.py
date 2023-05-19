import os
from multiprocessing import Manager, Pool

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
    clone, no_hla_data, hla_data, index = clone
    hla_clone = hla_data[clone].sum()
    N = get_sum_clones_by_runs(hla_data.run)
    P = (hla_data[clone].sum() + no_hla_data[clone].sum()) / get_sum_clones_by_runs()
    p_val = 1 - binom.cdf(hla_clone, N, P)
    clone_to_pval[clone] = p_val


def process_one_clone_fisher(clone):
    clone, no_hla_data, hla_data, index = clone
    hla_clone = hla_data[clone].sum()
    hla_no_clone = get_sum_clones_by_runs(hla_data.run) - hla_clone
    no_hla_clone = no_hla_data[clone].sum()
    no_hla_no_clone = get_sum_clones_by_runs(no_hla_data.run) - no_hla_clone
    res = fisher_exact([[hla_clone, hla_no_clone], [no_hla_clone, no_hla_no_clone]])
    clone_to_pval[clone] = res[1]
    if index % 50000 == 0:
        print(f'Processed {index} clones,', end='')


def hla_test(hla_data, no_hla_data, clonotype_matrix, fisher=True):
    arguments = [(x, no_hla_data, hla_data, i) for i, x in enumerate(clonotype_matrix.columns)]
    with Pool(40) as p:
        if fisher:
            p.map(process_one_clone_fisher, arguments)
        else:
            p.map(process_one_clone_binom, arguments)
    pvals = []
    for clone in tqdm(clonotype_matrix.columns):
        pvals.append(clone_to_pval[clone])
    significant_pvals = hochberg(pvals, alpha=0.05)
    significant_clones = []
    for pval, clone in zip(significant_pvals, clonotype_matrix.columns):
        if pval:
            significant_clones.append(clone)
    return significant_clones


def create_subdata_from_clonotype_matrix(clonotype_matrix, hla_desc):
    hla_data = clonotype_matrix[clonotype_matrix.run.isin(hla_desc)]
    no_hla_data = clonotype_matrix[~clonotype_matrix.run.isin(hla_desc)]
    clonotype_matrix = clonotype_matrix.drop(columns=['run'])
    return clonotype_matrix, hla_data, no_hla_data


def prepare_clonotype_matrix(clonotype_matrix_path, desc):
    return pd.read_csv(clonotype_matrix_path).drop_duplicates().set_index('cdr3aa').T.reset_index().rename(
        columns={'index': 'run'}).merge(desc[['run']])


def hla_test_allele_based(allele, hla_desc_path, clone_matrix_path, save_path, fisher=True):
    desc_fmba_not_nan = pd.read_csv(f'data/desc_fmba_not_nan_hla.csv').drop(columns=['Unnamed: 0'])
    hla_desc = set(pd.read_csv(f'{hla_desc_path}/fmba_desc_hla_{allele}.csv').run)

    clonotype_matrix = prepare_clonotype_matrix(
        f'{clone_matrix_path}/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{allele}.csv', desc_fmba_not_nan)
    print(clonotype_matrix.head())

    clonotype_matrix, hla_data, no_hla_data = create_subdata_from_clonotype_matrix(clonotype_matrix, hla_desc)
    print(hla_data.shape, no_hla_data.shape)
    significant_clones = hla_test(hla_data, no_hla_data, clonotype_matrix, fisher)
    pd.DataFrame(data={'clone': significant_clones}).to_csv(
        f'{save_path}/hla_associated_clones_500k_top_1_mismatch_hla_{allele}.csv',
        index=False)


if __name__ == "__main__":
    run_to_number_of_clonotypes = pd.read_csv(snakemake.input[0]).set_index('run')
    hla_keys = snakemake.params.hla_to_consider
    print(hla_keys)
    if not os.path.exists(snakemake.output[0]):
        os.mkdir(snakemake.output[0])
    for hla in list(hla_keys):
        hla_test_allele_based(hla,
                              hla_desc_path=snakemake.input[1],
                              clone_matrix_path=snakemake.input[2],
                              save_path=snakemake.output[0],
                              fisher=True)
