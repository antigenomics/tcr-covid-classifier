import random

import igraph as ig
import numpy as np
import pandas as pd
from multipy.fwer import hochberg
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact

random.seed(42)


def hdist(s1, s2):
    if len(s1) != len(s2):
        return float('inf')
    else:
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def seqs2hamming(seqs, threshold=1, viz_method='graphopt'):
    seqs = np.array(seqs).astype("str")
    dm = squareform(pdist(seqs.reshape(-1, 1), metric=lambda x, y: hdist(x[0], y[0])))
    dmf = pd.DataFrame(dm, index=seqs, columns=seqs).stack().reset_index()
    dmf.columns = ['id1', 'id2', 'distance']
    dmf = dmf[dmf['distance'] <= threshold]

    # graph
    graph = ig.Graph.TupleList(dmf[['id1', 'id2']].itertuples(index=False))

    # clusters
    clusters = graph.components()
    membership = clusters.membership

    # layout
    layout = graph.layout(viz_method)
    coords = np.array(layout.coords)

    df_graph = pd.DataFrame(
        {'cdr3': graph.vs()['name'],
         'cluster': membership,
         'x': coords[:, 0],
         'y': coords[:, 1]
         })

    # summary
    df_graph_summary = df_graph.groupby(['cluster']).agg(
        cluster_size=('cluster', 'size'),
        x_mean=('x', 'mean'),
        y_mean=('y', 'mean')).reset_index()

    return pd.merge(df_graph, df_graph_summary)


def get_network_components(seqs, threshold=1):
    seqs = np.array(seqs).astype("str")
    dm = squareform(pdist(seqs.reshape(-1, 1), metric=lambda x, y: hdist(x[0], y[0])))
    dmf = pd.DataFrame(dm, index=seqs, columns=seqs).stack().reset_index()
    dmf.columns = ['id1', 'id2', 'distance']
    dmf = dmf[dmf['distance'] <= threshold]

    # graph
    graph = ig.Graph.TupleList(dmf[['id1', 'id2']].itertuples(index=False))

    # clusters
    clusters = graph.components()
    return clusters


def check_distance(cdr1, cdr2, dist=1):
    if len(cdr1) != len(cdr2):
        return False
    found_diff_count = 0
    for c1, c2 in zip(cdr1, cdr2):
        if c1 != c2 and found_diff_count < dist:
            found_diff_count += 1
        elif c1 != c2 and found_diff_count == dist:
            return False
    return True


def check_db_epitopes_cdr3(db, cdr3, dist=1):
    return db[db.cdr3.apply(lambda x: check_distance(x, cdr3, dist))][['antigen.epitope', 'antigen.species']]


def get_epitopes_for_beta_clone(vdjdb, beta_cdr3, dist=1):
    beta_epitopes = check_db_epitopes_cdr3(vdjdb, beta_cdr3, dist=dist)
    return beta_epitopes.drop_duplicates()  # .merge(alpha_epitopes).drop_duplicates()


def get_epitopes_for_cluster(vdjdb, clones_to_cluster, cluster, dist=1):
    cur_cluster = clones_to_cluster[clones_to_cluster.cluster == cluster]
    all_data = []
    for cdr3 in cur_cluster.cdr3:
        all_data.append(get_epitopes_for_beta_clone(vdjdb, cdr3, dist=dist))
    res = pd.concat(all_data)
    res['count'] = 1
    return res.groupby(['antigen.epitope', 'antigen.species'], as_index=False).count().sort_values(by='count')


def get_count_of_antigen_associated_clones(vdjdb, antigen, chain='TRB'):
    return len(vdjdb[(vdjdb.gene == chain) & (vdjdb['antigen.epitope'] == antigen)]), len(vdjdb[(vdjdb.gene == chain)])


def check_significant_epitopes_for_cluster(vdjdb, res_beta, cluster, dist=1, gene='TRB', alpha=0.05):
    epitopes = get_epitopes_for_cluster(vdjdb, res_beta, cluster, dist)
    trb_in_vdjdb = len(vdjdb[vdjdb.gene.str.contains(gene)]['antigen.epitope'])
    pvals = []
    if len(epitopes) == 0:
        return None
    for epi, count in zip(epitopes['antigen.epitope'], epitopes['count']):
        x = get_count_of_antigen_associated_clones(vdjdb, epi, chain=gene)[0]
        y = count
        pvals.append(fisher_exact([[x, trb_in_vdjdb - x], [y, len(res_beta) - y]], alternative='less')[1])
    if len(pvals) > 1:
        sign = hochberg(pvals, alpha=alpha)
    else:
        sign = [pvals[0] < alpha]
    return epitopes[sign] if len(epitopes[sign]) > 0 else None


def get_epitope_for_clone(cdr3, vdjdb):
    res = get_epitopes_for_beta_clone(vdjdb, cdr3, dist=1)
    if len(res) > 0:
        return list(res['antigen.epitope'])[0]
    return None


def get_significant_epitopes_to_clone_mapping(vdjdb, res_beta, cluster, significant_epitopes_for_cluster, dist=1,
                                              gene='TRB'):
    result = res_beta[res_beta.cluster == cluster]
    gene_vdjdb = vdjdb[
        (vdjdb.gene.str.contains(gene)) & (vdjdb['antigen.epitope'].isin(significant_epitopes_for_cluster))]
    result['epitope'] = result.cdr3.apply(lambda x: get_epitope_for_clone(x, gene_vdjdb))
    return result


def get_most_frequent_cluster_by_vdjdb_occurence(vdjdb, cluster_epitopes, gene='TRB'):
    cluster_epitopes['cluster_epitopes_freq'] = cluster_epitopes.apply(lambda x: x['count'] / vdjdb[
        (vdjdb.gene == gene) & (vdjdb['antigen.epitope'] == x['antigen.epitope'])].cdr3.nunique(), axis=1)
    return cluster_epitopes.sort_values(by='cluster_epitopes_freq', ascending=False).loc[0, :]


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


def plot_clusters_of_clonotypes(clustering_res, color_by='cluster', ax=None):
    if ax is None:
        fig, (ax) = plt.subplots(1, 1)
    if len(clustering_res[color_by].unique()) > 10:
        palette = sns.color_palette("tab20b", 100)
    else:
        palette = sns.color_palette("tab10")
    sns.scatterplot(clustering_res[clustering_res.cluster_size > 1], x='x', y='y', hue=color_by,
                    palette=palette, ax=ax)
    sns.scatterplot(clustering_res[clustering_res.cluster_size == 1], x='x', y='y', hue=color_by, palette=['grey'],
                    legend=False, ax=ax)
