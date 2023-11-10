import os,sys,inspect
print(os.getcwd())
sys.path.append(os.getcwd())

import pandas as pd
from utils.clustering_utils import seqs2hamming
from utils.data_utils import prepare_clonotype_matrix


def identify_cluster_presence_in_sample(cm):
    return cm.sum(axis=1).astype(bool).astype(int)


if __name__ == "__main__":
    gene = snakemake.params.gene
    cm_path = snakemake.input[0]
    matrix = prepare_clonotype_matrix(cm_path)
    clustering_res = seqs2hamming(matrix.columns[1:], viz_method='drl')
    cluster_matrix = pd.DataFrame(data={'run': matrix['run']})
    for cluster in clustering_res.cluster.unique():
        cluster_matrix[f'{gene}_cluster_{cluster}'] = identify_cluster_presence_in_sample(
            matrix[clustering_res[clustering_res.cluster == cluster].cdr3])
    cluster_matrix.to_csv(snakemake.output[0])
