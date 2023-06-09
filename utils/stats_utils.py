import scipy.stats as stats
from scipy.stats import chi2_contingency
from utils.data_utils import prepare_clonotype_matrix, prepare_run_column
import pandas as pd
import numpy as np


def evaluate_anova_testing(matrix, by, pvalue_threshold=0.05):
    important_list = []
    names = list(matrix[by].unique())
    v_genes = list(x for x in matrix.columns if x.startswith('TR'))
    for v_gene in v_genes:
        datasets = [matrix[v_gene][matrix[by] == name] for name in names]
        pvalue = stats.f_oneway(*datasets)[1]
        if pvalue < pvalue_threshold:
            important_list.append(v_gene)
            print(f'There is a significant difference for {v_gene} broken by {by} (p-value={pvalue})')
    return important_list


def evaluate_mannwhitneyu_testing(matrix, by, pvalue_threshold=0.05):
    important_list = []
    names = list(matrix[by].unique())
    v_genes = list(x for x in matrix.columns if x.startswith('TR'))
    for v_gene in v_genes:
        datasets = [matrix[v_gene][matrix[by] == name] for name in names]
        pvalue = stats.mannwhitneyu(*datasets)[1]
        if pvalue < pvalue_threshold:
            important_list.append(v_gene)
            print(f'There is a significant difference for {v_gene} broken by {by} (p-value={pvalue})')
    return important_list


def get_sum_clones_by_runs(run_to_number_of_clonotypes, runs=None):
    if runs is None:
        return run_to_number_of_clonotypes.number_of_clones.sum()
    return run_to_number_of_clonotypes.loc[runs, :].number_of_clones.sum()


def evaluate_fisher(clone, no_feature_data, feature_data, run_to_number_of_clonotypes):
    feature_clone = feature_data[clone].sum()
    feature_no_clone = get_sum_clones_by_runs(run_to_number_of_clonotypes, feature_data.run) - feature_clone
    no_feature_clone = no_feature_data[clone].sum()
    no_feature_no_clone = get_sum_clones_by_runs(run_to_number_of_clonotypes, no_feature_data.run) - no_feature_clone
    res = chi2_contingency([[feature_clone, feature_no_clone], [no_feature_clone, no_feature_no_clone]])
    if sum(map(lambda x: x < 5, res[3].flatten())) != 0:
        res = stats.fisher_exact([[feature_clone, feature_no_clone], [no_feature_clone, no_feature_no_clone]])
    return res[1]


def get_top_changed_clonotypes(clonotype_matrix, desc, pvals, log_fold_change_threshold, logp_threshold,
                               healthy_col='covid', healthy_label='healthy'):
    print('v')
    cm = prepare_run_column(clonotype_matrix).merge(
        prepare_run_column(desc[['run', healthy_col]])
    )

    healthy_data = cm[cm[healthy_col] == healthy_label].drop_duplicates().set_index('run').drop(
        columns=[healthy_col]).T.reset_index().rename(
        columns={'index': 'clone'}).set_index('clone')
    healthy_data['count_of_ways_h'] = healthy_data.sum(axis=1)
    healthy_data = healthy_data.reset_index()[['count_of_ways_h', 'clone']]

    ill_data = cm[cm[healthy_col] != healthy_label].drop_duplicates().set_index('run').drop(
        columns=[healthy_col]).T.reset_index().rename(
        columns={'index': 'clone'}).set_index('clone')
    ill_data['count_of_ways_i'] = ill_data.sum(axis=1)
    ill_data = ill_data.reset_index()[['count_of_ways_i', 'clone']]

    df = pvals.merge(healthy_data).merge(ill_data)
    df['fold_change'] = df['count_of_ways_i'] / df['count_of_ways_h']
    df['log_fold_change'] = df['fold_change'].apply(lambda x: np.log2(x))
    return df