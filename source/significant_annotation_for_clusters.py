import os
import sys

import numpy as np
import pandas as pd
from multipy.fdr import lsu

print(os.getcwd())
sys.path.append(os.getcwd())

from utils.clustering_utils import seqs2hamming, check_significant_epitopes_for_all_clusters


def create_significant_associations_df(cm_path, vdjdb_path, save_path, prefix='fmba', gene='TRB'):
    covid_clones_beta = []
    for cm in cm_path:
        clean_beta_cm = pd.read_csv(cm).drop(columns=['Unnamed: 0'])
        covid_clones_beta += list(clean_beta_cm.columns[1:])
    res_beta = seqs2hamming(covid_clones_beta, viz_method='drl')
    vdjdb = pd.read_csv(vdjdb_path, sep='\t')
    all_associations, epitopes = check_significant_epitopes_for_all_clusters(res=res_beta,
                                                                   vdjdb=vdjdb,
                                                                   gene=gene,
                                                                   alpha=0.05,
                                                                   threads=snakemake.threads,
                                                                   dir_to_save=None)
    associations_per_cluster = []
    for k in all_associations:
        cur_cluster_associations = all_associations[k]
        if cur_cluster_associations is not None:
            cur_cluster_associations['cluster'] = k
            associations_per_cluster.append(cur_cluster_associations)
    pvals_result_beta = pd.concat(associations_per_cluster)
    pvals_result_beta = pvals_result_beta[lsu(np.array(pvals_result_beta.pval), q=0.01)]
    pvals_result_beta['all_clust_count'] = pvals_result_beta['antigen.epitope'].apply(lambda x: epitopes[x])
    pvals_result_beta['enrichment_score'] = pvals_result_beta['count'] / pvals_result_beta['all_clust_count']
    pvals_result_beta = pvals_result_beta[pvals_result_beta['antigen.species'] != 'HomoSapiens']
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    pvals_result_beta.to_csv(f'{save_path}/{prefix}_{gene}_sign_assoc.csv')

    max_enrichment_df = []
    for cluster in sorted(pvals_result_beta.cluster.unique()):
        max_enrichment_score = pvals_result_beta[pvals_result_beta.cluster == cluster].enrichment_score.max()
        max_enrichment_df.append(pvals_result_beta[(pvals_result_beta.cluster == cluster) & (
                pvals_result_beta.enrichment_score == max_enrichment_score)])
    pd.concat(max_enrichment_df).to_csv(f'{save_path}/{prefix}_{gene}_sign_assoc_with_max_enrichment.csv')


if __name__ == "__main__":
    if 'snakemake' in globals():
        create_significant_associations_df(cm_path=snakemake.input[0:-1],
                                           vdjdb_path=snakemake.input[-1],
                                           save_path=snakemake.params.working_directory,
                                           prefix=snakemake.params.prefix,
                                           gene=snakemake.params.gene)
