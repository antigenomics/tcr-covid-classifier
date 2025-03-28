import os
import random
from multiprocessing import Manager, Pool

import numpy as np
import pandas as pd
from multipy.fdr import lsu
from multipy.fwer import hochberg
from scipy.stats import fisher_exact, binom
from tqdm import tqdm

clone_to_pval = Manager().dict()


def prepare_run_column(df):
    df['run'] = df['run'].apply(lambda x: x.split('.')[0])
    return df


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


def process_one_clone_fisher(clone):
    clone, healthy_data, covid_data, alternative = clone
    covid_clone = covid_data[clone].sum()
    covid_no_clone = get_sum_clones_by_runs(covid_data.run) - covid_clone
    no_covid_clone = healthy_data[clone].sum()
    no_covid_no_clone = get_sum_clones_by_runs(healthy_data.run) - no_covid_clone
    # res = chi2_contingency([[covid_clone, covid_no_clone], [no_covid_clone, no_covid_no_clone]])
    # if sum(map(lambda x: x < 5, res[3].flatten())) != 0:
    res = fisher_exact([[covid_clone, covid_no_clone], [no_covid_clone, no_covid_no_clone]], alternative=alternative)

    clone_to_pval[clone] = res[1]


def covid_test(healthy_data, covid_data, clonotype_matrix, pval_save_path=None, fisher=True, alternative='greater'):
    arguments = [(x, healthy_data[['run', x]], covid_data[['run', x]], alternative) for x in clonotype_matrix.columns]
    print(f'Testing started')
    with Pool(snakemake.threads, maxtasksperchild=2) as p:
        if fisher:
            p.map(process_one_clone_fisher, arguments)
            print('finished')
        else:
            p.map(process_one_clone_binom, arguments)
    pvals = []
    for clone in tqdm(clonotype_matrix.columns):
        pvals.append(clone_to_pval[clone])
    if 'fmba' in snakemake.params.platform:
        significant_pvals = lsu(np.array(pvals), q=significant_threshold)
        for q in [0.005, 0.01, 0.05, 0.1]:
            print(f'There are {sum(lsu(np.array(pvals), q=q))} for threshold={q}')
    else:
        significant_pvals = hochberg(pvals, alpha=significant_threshold)
    if pval_save_path is not None:
        pd.DataFrame(data={'clone': clonotype_matrix.columns, 'pval': pvals}).to_csv(pval_save_path)
    significant_clones = []
    for pval, clone in zip(significant_pvals, clonotype_matrix.columns):
        if pval:
            significant_clones.append(clone)
    return significant_clones


def covid_test_matrix_based(clonotype_matrix, desc_path, save_path,
                            pval_save_path=None, n=None, platform=None,
                            fisher=True, marker_column='covid',
                            marker_column_success_sign='covid', alternative='greater', drop_test=False):
    desc = pd.read_csv(desc_path).drop(columns=['Unnamed: 0'])
    if 'is_test_run' in desc.columns and drop_test:
        print('test was dropped')
        desc = desc[~desc.is_test_run]
    if platform is not None:
        desc = desc[desc.platform == platform]
    if n is not None:
        clonotype_matrix = clonotype_matrix.head(n)
    if 'cdr3aa' in clonotype_matrix.columns:
        clonotype_matrix = clonotype_matrix.drop_duplicates().set_index('cdr3aa').T.reset_index().rename(
            columns={'index': 'run'})
    if 'Unnamed: 0' in clonotype_matrix.columns:
        clonotype_matrix = clonotype_matrix.drop(columns=['Unnamed: 0'])
    clonotype_matrix = prepare_run_column(clonotype_matrix).merge(prepare_run_column(desc[['run', marker_column]]))
    if sample_random is not None:
        random.seed(42)
        sampled_runs = random.sample(set(clonotype_matrix.run), sample_random)
        clonotype_matrix = clonotype_matrix[clonotype_matrix.run.isin(sampled_runs)]
    print(clonotype_matrix.shape)
    covid_data = clonotype_matrix[clonotype_matrix[marker_column] == marker_column_success_sign]
    healthy_data = clonotype_matrix[clonotype_matrix[marker_column] != marker_column_success_sign]
    clonotype_matrix = clonotype_matrix.drop(columns=[marker_column, 'run'])
    print(f'Finished preparing data. Healthy data shape: {healthy_data.shape}, ill: {covid_data.shape}')

    significant_clones = covid_test(healthy_data, covid_data, clonotype_matrix, pval_save_path, fisher, alternative)

    pd.DataFrame(data={'clone': significant_clones}).to_csv(save_path, index=False)


def covid_test_allele_based(save_path, cm, desc_path, pval_save_path, fisher, drop_test, alternative='two-sided'):
    covid_test_matrix_based(clonotype_matrix=cm,
                            desc_path=desc_path,
                            save_path=save_path,
                            pval_save_path=pval_save_path,
                            fisher=fisher,
                            drop_test=drop_test,
                            alternative=alternative)


def covid_test_for_beta_chain_adaptive(n=500000):
    clm = pd.read_csv('data/adaptive_clone_results/clonotype_matrix_fmba_adaptive_top_500k_1_mismatch.csv')
    print(f'Finished reading clonotype matrix')
    covid_test_matrix_based(
        clonotype_matrix=clm,
        desc_path=f'data/standardized_log_exp_usage_matrix_joint_new.csv',
        save_path=f'data/adaptive_clone_results/covid_clones_all_hla_top_500k_1_mismatch_fisher.csv',
        pval_save_path=f'data/adaptive_clone_results/covid_pvals_all_hla_top_500k_1_mismatch_fisher.csv',
        n=n,
        fisher=True,
        platform='adaptive'
    )


def batch_test_for_beta_chain_adaptive(n=500000):
    clm = pd.read_csv('data/adaptive_clone_results/clonotype_matrix_fmba_adaptive_top_500k_1_mismatch.csv')
    print(f'Finished reading clonotype matrix')
    covid_test_matrix_based(
        clonotype_matrix=clm,
        desc_path=f'data/desc_hip_bool.csv',
        save_path=f'data/adaptive_clone_results/batch_clones_top_500k_1_mismatch_fisher.csv',
        pval_save_path=f'data/adaptive_clone_results/batch_pvals_top_500k_1_mismatch_fisher.csv',
        n=n,
        fisher=True,
        platform='adaptive',
        marker_column='is_keck',
        marker_column_success_sign='yes',
        alternative='two-sided'
    )


def covid_test_for_alpha_chain(n=500000):
    global run_to_number_of_clonotypes
    run_to_number_of_clonotypes = pd.read_csv('data/run_to_number_of_clonotypes_alpha.csv').set_index('run')

    covid_test_matrix_based(
        clonotype_matrix=pd.read_csv('data/alpha/clonotype_matrix_fmba_alpha_top_500k_1_mismatch.csv'),
        desc_path=f'data/standardized_log_exp_usage_matrix_by_v_gene_alpha_fmba.csv',
        save_path=f'data/alpha/covid_clones_all_hla_top_{n // 1000}k_1_mismatch_fisher.csv',
        pval_save_path=f'data/alpha/clone_pvals__all_hla_top_{n // 1000}k_1_mismatch_fisher.csv',
        n=n,
        fisher=True)


def covid_test_for_fmba(run_to_number_of_clones_path, clonotype_matrix_path, um_path, save_path, pval_save_path,
                        drop_test):
    global run_to_number_of_clonotypes
    run_to_number_of_clonotypes = prepare_run_column(pd.read_csv(run_to_number_of_clones_path)).set_index('run')

    covid_test_matrix_based(
        clonotype_matrix=pd.read_csv(clonotype_matrix_path),
        desc_path=um_path,
        save_path=save_path,
        pval_save_path=pval_save_path,
        fisher=True,
        drop_test=drop_test
    )


def test_for_all_hla_alleles(run_to_number_of_clones_path, cm_folder_path, desc_folder_path, save_path, fisher,
                             drop_test, alternative, hla_keys=None):
    global run_to_number_of_clonotypes
    run_to_number_of_clonotypes = prepare_run_column(pd.read_csv(run_to_number_of_clones_path)).set_index('run')

    if hla_keys is None:
        hla_keys = pd.read_csv('data/hla_keys.csv')['0']
    for hla in list(hla_keys):
        print(f'started processing {hla}')
        cm = pd.read_csv(f'{cm_folder_path}/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{hla}.csv')
        desc_path = f'{desc_folder_path}/fmba_desc_hla_{hla}.csv'
        covid_test_allele_based(save_path + f'/hla_covid_associated_clones_500k_top_1_mismatch_hla_{hla}.csv',
                                cm,
                                desc_path,
                                save_path + f'/hla_covid_associated_clones_1_mismatch_hla_{hla}_pvals.csv',
                                fisher=fisher,
                                drop_test=drop_test,
                                alternative=alternative
                                )


if __name__ == "__main__":
    if 'snakemake' in globals():
        significant_threshold = snakemake.params.significant_threshold
        print(f'Multiple test threshold {significant_threshold}')
        drop_test = snakemake.params.drop_test
        print(f'Test to be dropped: {drop_test}')
        platform = snakemake.params.platform.split('-')[0]
        sample_random = None
        if 'random' in snakemake.params.platform.split('-'):
            sample_random = int(snakemake.params.platform.split('-')[-1])
        if platform == 'fmba' or platform == 'adaptive':
            covid_test_for_fmba(run_to_number_of_clones_path=snakemake.input[1],
                                clonotype_matrix_path=snakemake.input[2],
                                um_path=snakemake.input[0],
                                save_path=snakemake.output[0],
                                pval_save_path=snakemake.output[1],
                                drop_test=drop_test
                                )
        if platform == 'vdjdb':
            run_to_number_of_clonotypes = prepare_run_column(
                pd.read_csv(snakemake.input[1])).set_index('run')

            covid_test_matrix_based(
                clonotype_matrix=pd.read_csv(snakemake.input[2]),
                desc_path=snakemake.input[0],
                save_path=snakemake.output[0],
                pval_save_path=snakemake.output[1],
                fisher=True,
                drop_test=drop_test,
                alternative='two-sided'
            )
        if platform == 'allele':
            print('making allele testing')
            if snakemake.params.hla_to_consider != []:
                hla_keys = snakemake.params.hla_to_consider
            else:
                hla_keys = pd.read_csv('data/hla_keys.csv')['0']
            print(hla_keys)
            if not os.path.exists(snakemake.output[0]):
                os.mkdir(snakemake.output[0])
            test_for_all_hla_alleles(hla_keys=hla_keys,
                                     run_to_number_of_clones_path=snakemake.input[0],
                                     cm_folder_path=snakemake.input[2],
                                     desc_folder_path=snakemake.input[1],
                                     save_path=snakemake.output[0],
                                     fisher=True,
                                     drop_test=drop_test,
                                     alternative='two-sided')
