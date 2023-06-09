import os

from utils.clustering_utils import check_significant_epitopes_for_cluster, seqs2hamming

import numpy as np
import pandas as pd
from tqdm import tqdm
from multiprocessing import Manager, Pool
from utils.data_utils import prepare_clonotype_matrix
from utils.viz_utils import plot_pandas_df_into_png
clone_to_df = Manager().dict()


def save_current_data(save_path):
    datasets_to_concat = []
    for clone, found_data in clone_to_df.items():
        datasets_to_concat.append(found_data)

    pd.concat(datasets_to_concat).to_csv(save_path)


def find_associated_alpha_clones_for_beta_clone(beta_clone, beta_cm, alpha_cm, save_path, alpha=0.95):
    runs_with_beta_clone = beta_cm[beta_cm[beta_clone] != 0].run
    runs_without_beta_clone = beta_cm[beta_cm[beta_clone] == 0].run
    alpha_cm_for_beta_clone = alpha_cm[alpha_cm.run.isin(runs_with_beta_clone)].drop(columns=['run'])
    usages = alpha_cm_for_beta_clone.sum(axis=0) / len(runs_with_beta_clone)
    usages = pd.DataFrame(usages, columns=['usage']).reset_index().rename(columns={'index': 'alpha_clone'})
    usages = usages[usages.usage > alpha]
    usages['beta_clone'] = beta_clone
    clone_to_df[beta_clone] = usages
    clone_to_df_size = len(clone_to_df)
    print(f'Processed {beta_clone}, iteration {clone_to_df_size}, found {len(usages)} associated alphas')
    if clone_to_df_size % 1000 == 0:
        save_current_data(save_path)
    return usages


def mapping_function(args):
    beta_clone, beta_cm, alpha_cm, save_path, alpha = args
    return find_associated_alpha_clones_for_beta_clone(beta_clone, beta_cm, alpha_cm, alpha)


def create_pairs_data(alpha_clone_matrix_path, beta_clone_matrix_path, save_path, alpha=0.95):
    datasets_to_concat = []
    alpha_cm = prepare_clonotype_matrix(alpha_clone_matrix_path, make_bool_features=True)
    beta_cm = prepare_clonotype_matrix(beta_clone_matrix_path, make_bool_features=True)
    print(f'Started searching for associations')
    arguments = []
    for i, clone in tqdm(enumerate(beta_cm.columns)):
        if clone == 'run':
            continue
        arguments.append([clone, beta_cm, alpha_cm, save_path, alpha])

    with Pool(40) as p:
        p.map(mapping_function, arguments)

    for clone, found_data in clone_to_df.items():
        datasets_to_concat.append(found_data)

    pd.concat(datasets_to_concat).to_csv(save_path)


def get_cooccurence_value_for_clusters(alpha_cdrs, beta_cdrs, alpha_matrix, beta_matrix, pairing_param=0.8):
    all_counter = 0
    success_counter = 0
    for alpha_clone in alpha_cdrs:
        for beta_clone in beta_cdrs:
            cur_df = pd.DataFrame({'alpha': alpha_matrix[alpha_clone], 'beta': beta_matrix[beta_clone]})
            cur_df['together'] = cur_df.alpha.astype(bool) & cur_df.beta.astype(bool)
            if sum(cur_df['together']) / cur_df.shape[0] > pairing_param:
                success_counter += 1
                break
        all_counter += 1
    return success_counter / all_counter


def make_metaclone_cm(cm, cluster_info):
    new_cm_cols = {'run': cm['run']}
    for cluster in cluster_info.cluster.unique():
        clones = cluster_info[cluster_info.cluster == cluster].cdr3
        new_cm_cols[f'cluster_{cluster}'] = cm[clones].sum(axis=1)
    return pd.DataFrame(new_cm_cols)


def plot_cooccured_epitopes_table(res_alpha, res_beta, cooccurence_dist1_epitopes, vdjdb, save_path=None, ax=None):
    alpha_index = []
    beta_index = []
    epi = []
    species = []
    for beta_cluster in res_beta.cluster.unique():
        for alpha_cluster in res_alpha.cluster.unique():
            if len(cooccurence_dist1_epitopes[beta_cluster][alpha_cluster]) > 0:
                for e in cooccurence_dist1_epitopes[beta_cluster][alpha_cluster]:
                    alpha_index.append(alpha_cluster)
                    beta_index.append(beta_cluster)
                    epi.append(e)
                    species.append(list(vdjdb[vdjdb['antigen.epitope'] == e]['antigen.species'])[0])

    df = pd.DataFrame(data={'alpha_cluster': alpha_index, 'beta_cluster': beta_index,
                            'epitope': epi, 'antigen_species': species})
    plot_pandas_df_into_png(df, output_path=save_path, ax=ax)


def evaluate_value_for_alpha_beta(alpha, beta):
    beta_ind = f'TRB_cluster_{beta}'
    alpha_ind = f'TRA_cluster_{alpha}'
    df = beta_cluster_presence[['run', beta_ind]].merge(
            alpha_cluster_presence[['run', alpha_ind]]).drop(
                columns=['run'])
    df['joint'] = df.apply(lambda x: 1 if x[beta_ind] + x[alpha_ind] == 2 else 0, axis=1)
    return df['joint'].sum() / (df[beta_ind].sum() + df[alpha_ind].sum() - df['joint'].sum())


def alpha_beta_joint_usage_matrix_preparation(tra_cm_path, trb_cm_path, vdjdb_path):
    vdjdb = pd.read_csv(vdjdb_path, sep='\t')
    alpha_matrix = pd.read_csv(tra_cm_path).drop(columns=['Unnamed: 0'])
    beta_matrix = pd.read_csv(trb_cm_path).drop(columns=['Unnamed: 0'])

    res_alpha = seqs2hamming(alpha_matrix.columns[1:], viz_method='drl')
    res_beta = seqs2hamming(beta_matrix.columns[1:], viz_method='graphopt')

    cooccurence90 = [[0 for _ in range(res_alpha.cluster.max() + 1)] for _ in range(res_beta.cluster.max() + 1)]
    for alpha_cluster in tqdm(range(res_alpha.cluster.max() + 1)):
        for beta_cluster in range(res_beta.cluster.max() + 1):
            cooccurence90[beta_cluster][alpha_cluster] = evaluate_value_for_alpha_beta(alpha_cluster, beta_cluster)

    beta_epitopes_dist_1 = {}
    for cluster in tqdm(res_beta.cluster.unique()):
        test_results = check_significant_epitopes_for_cluster(vdjdb, res_beta, cluster, dist=1, gene='TRB')
        beta_epitopes_dist_1[cluster] = set(test_results['antigen.epitope']) if test_results is not None else set()

    alpha_epitopes_dist_1 = {}
    for cluster in tqdm(res_alpha.cluster.unique()):
        test_results = check_significant_epitopes_for_cluster(vdjdb, res_alpha, cluster, dist=1, gene='TRA')
        alpha_epitopes_dist_1[cluster] = set(test_results['antigen.epitope']) if test_results is not None else set()

    cooccurence_dist1_epitopes = [[set() for _ in range(res_alpha.cluster.max() + 1)] for _ in
                                  range(res_beta.cluster.max() + 1)]
    for beta_cluster in res_beta.cluster.unique():
        for alpha_cluster in res_alpha.cluster.unique():
            cooccurence_dist1_epitopes[beta_cluster][alpha_cluster] = alpha_epitopes_dist_1[alpha_cluster].intersection(
                beta_epitopes_dist_1[beta_cluster])

    plotting_df = pd.DataFrame(data={i: x for i, x in enumerate(np.array(cooccurence90))})
    annot_df = pd.DataFrame(data=np.array([[len(x) for x in y] for y in cooccurence_dist1_epitopes]).T)

    plotting_df.to_csv(snakemake.output[0],
                       index=False)
    annot_df.to_csv(snakemake.output[1], index=False)

    make_metaclone_cm(alpha_matrix, res_alpha).to_csv(snakemake.output[2])
    make_metaclone_cm(beta_matrix, res_beta).to_csv(snakemake.output[3])
    plot_cooccured_epitopes_table(res_alpha, res_beta, cooccurence_dist1_epitopes, vdjdb, snakemake.output[4])


if __name__ == "__main__":
    beta_cluster_presence = pd.read_csv(snakemake.input[3]).drop(columns=['Unnamed: 0'])
    alpha_cluster_presence = pd.read_csv(snakemake.input[4]).drop(columns=['Unnamed: 0'])
    alpha_beta_joint_usage_matrix_preparation(tra_cm_path=snakemake.input[0],
                                              trb_cm_path=snakemake.input[1],
                                              vdjdb_path=snakemake.input[2],
                                              )
    # create_pairs_data(alpha_clone_matrix_path='data/alpha/clonotype_matrix_fmba_alpha_top_500k_0_mismatch.csv',
    #                   beta_clone_matrix_path='data/hla_sign_clone_matrix/covid_clonotype_matrix_500k_top_1_mismatch_hla_all_fisher.csv',
    #                   save_path='data/alpha_beta_pairs/paired_clones_alpha_mis_0_beta_mis_1_alpha_095_covid.csv',
    #                   alpha=0.95)
