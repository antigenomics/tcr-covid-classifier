import random

import igraph as ig
import numpy as np
import pandas as pd
from multipy.fwer import hochberg
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact
from multipy.fdr import lsu
from multiprocessing import Pool, Manager
from tqdm import trange

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
    return db[db.cdr3.apply(lambda x: check_distance(x, cdr3, dist))][['antigen.epitope', 'antigen.species']].drop_duplicates()


def get_epitopes_for_beta_clone(vdjdb, beta_cdr3, dist=1):
    # print(len(vdjdb))
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


def get_count_of_antigen_associated_clones(vdjdb, antigen, chain='TRB', res_beta=None):
    number_of_matches = 0
    if res_beta is not None:
        for clone in res_beta.cdr3:
            if antigen in set(get_epitopes_for_beta_clone(vdjdb[vdjdb['antigen.epitope'] == antigen], clone)['antigen.epitope']):
                number_of_matches += 1
        return number_of_matches
    return len(vdjdb[(vdjdb.gene == chain) & (vdjdb['antigen.epitope'] == antigen)]), len(vdjdb[(vdjdb.gene == chain)])


def check_significant_epitopes_for_cluster(vdjdb, res_beta, cluster, known_epitopes, dist=1, gene='TRB', alpha=0.05, threads=32):
    global check_one_epitope_significance
    # TODO
    def check_one_epitope_significance(args):
        epi, count, overall_trb, cluster_trb, gene= args
        if epi in known_epitopes:
            x = known_epitopes[epi]
        else:
            x = get_count_of_antigen_associated_clones(vdjdb, epi, res_beta=res_beta, chain=gene)
            known_epitopes[epi] = x
        y = count
        # print([[x, overall_trb - x], [y, cluster_trb - y]])
        assert y <= x
        pvals[epi] = fisher_exact([[x, overall_trb - x], [y, cluster_trb - y]], alternative='less')[1]

    epitopes = get_epitopes_for_cluster(vdjdb, res_beta, cluster, dist)
    overall_trb = len(res_beta)
    cluster_trb = len(res_beta[res_beta.cluster == cluster])
    pvals = Manager().dict()
    if len(epitopes) == 0:
        return None
    arguments = []
    for epi, count in zip(epitopes['antigen.epitope'], epitopes['count']):
        arguments.append([epi, count, overall_trb, cluster_trb, gene])
    with Pool(threads, maxtasksperchild=2) as p:
        p.map(check_one_epitope_significance, arguments)
    pvals_res = []
    for epi in epitopes['antigen.epitope']:
        pvals_res.append(pvals[epi])
    epitopes['pval'] = pd.Series(pvals_res)
    if len(pvals_res) > 1:
        sign = lsu(np.array(pvals_res), q=alpha)
    else:
        sign = [pvals_res[0] < alpha]
    return epitopes[sign] if len(epitopes[sign]) > 0 else None


def check_significant_epitopes_for_all_clusters(res, vdjdb, gene, alpha, threads, dir_to_save='data/fmba_associations'):
    # print('newnew!')
    cluster_to_epi = {}
    known_epitopes = Manager().dict()
    for cluster_index in trange(res.cluster.max() + 1):
        cluster_to_epi[cluster_index] = check_significant_epitopes_for_cluster(vdjdb, res, cluster_index,
                                                                                    dist=1, gene=gene, alpha=alpha,
                                                                                    threads=threads,
                                                                                    known_epitopes=known_epitopes)
        if cluster_to_epi[cluster_index] is not None and dir_to_save is not None:
            cluster_to_epi[cluster_index].to_csv(f'{dir_to_save}/{gene}_cluster_{cluster_index}.csv')
    if dir_to_save is not None:
        pd.DataFrame.from_dict(known_epitopes, orient='index').reset_index().to_csv(f'{dir_to_save}/{gene}_epitopes.csv')
    return cluster_to_epi, known_epitopes


def check_significant_epitopes_for_cluster_vdjdb_based(vdjdb, res_beta, cluster, dist=1, gene='TRB', alpha=0.05):
    # TODO
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
        sign = lsu(np.array(pvals), q=alpha)
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


def cluster_epitopes_freq_calc(x, vdjdb, gene):
    vdjdb_value = vdjdb[
        (vdjdb.gene == gene) & (vdjdb['antigen.epitope'] == x['antigen.epitope'])].cdr3.nunique()
    return 0 if vdjdb_value == 0 else x['count'] / vdjdb_value


def get_most_frequent_cluster_by_vdjdb_occurence(vdjdb, cluster_epitopes, gene='TRB'):
    cluster_epitopes['cluster_epitopes_freq'] = cluster_epitopes.apply(
        lambda x: cluster_epitopes_freq_calc(x, vdjdb, gene), axis=1)
    return cluster_epitopes.sort_values(by='cluster_epitopes_freq', ascending=False).reset_index(drop=True).loc[0, :]


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


def create_summary_stats_table(clustering_res, cluster_to_epi, cm, vdjdb, test_batch, desc, gene='TRB'):
    summary = clustering_res[['cluster', 'cluster_size']].drop_duplicates()
    summary['antigen.epitope'] = summary.cluster.apply(lambda x: get_most_frequent_cluster_by_vdjdb_occurence(
        vdjdb,
        cluster_to_epi[x][~cluster_to_epi[x]['antigen.species'].str.contains('apiens')].reset_index(drop=True),
        gene='TRB')[['antigen.epitope', 'antigen.species']][0]
                                                       )
    summary['antigen.species'] = summary.cluster.apply(lambda x: get_most_frequent_cluster_by_vdjdb_occurence(
        vdjdb,
        cluster_to_epi[x][~cluster_to_epi[x]['antigen.species'].str.contains('apiens')].reset_index(drop=True),
        gene=gene)[['antigen.epitope', 'antigen.species']][0]
                                                       )
    summary['epitope.num_clones'] = summary.apply(
        lambda x:
        cluster_to_epi[x.cluster][cluster_to_epi[x.cluster]['antigen.epitope'] == x['antigen.epitope']]['count'][0],
        axis=1)

    train_runs = desc[desc.project.str.contains(test_batch)]
    test_runs = desc[~desc.project.str.contains(test_batch)]

    healthy_train_cm = cm[cm.run.isin(train_runs)].merge(desc[['run']][desc.COVID_status == 'healthy'])
    healthy_test_cm = cm[cm.run.isin(test_runs)].merge(desc[['run']][desc.COVID_status == 'healthy'])
    covid_train_cm = cm[cm.run.isin(train_runs)].merge(desc[['run']][desc.COVID_status != 'healthy'])
    covid_test_cm = cm[cm.run.isin(test_runs)].merge(desc[['run']][desc.COVID_status != 'healthy'])

    train_cm = cm[cm.run.isin(train_runs)]
    test_cm = cm[cm.run.isin(test_runs)]

    summary['num_samples_with_cluster_train'] = summary.cluster.apply(
        lambda x: (train_cm[clustering_res[clustering_res.cluster == x].cdr3].sum(axis=1) > 0).sum())
    summary['num_samples_with_cluster_train_healthy'] = summary.cluster.apply(
        lambda x: (healthy_train_cm[clustering_res[clustering_res.cluster == x].cdr3].sum(axis=1) > 0).sum())
    summary['num_samples_with_cluster_train_covid'] = summary.cluster.apply(
        lambda x: (covid_train_cm[clustering_res[clustering_res.cluster == x].cdr3].sum(axis=1) > 0).sum())

    summary['num_samples_with_cluster_test'] = summary.cluster.apply(
        lambda x: (test_cm[clustering_res[clustering_res.cluster == x].cdr3].sum(axis=1) > 0).sum())
    summary['num_samples_with_cluster_test_healthy'] = summary.cluster.apply(
        lambda x: (healthy_test_cm[clustering_res[clustering_res.cluster == x].cdr3].sum(axis=1) > 0).sum())
    summary['num_samples_with_cluster_test_covid'] = summary.cluster.apply(
        lambda x: (covid_test_cm[clustering_res[clustering_res.cluster == x].cdr3].sum(axis=1) > 0).sum())

    summary.to_excel(f'figures/clustering_summary_{gene}.xlsx')
    summary.to_csv(f'figures/clustering_summary_{gene}.csv')


def read_association_data(path):
    df = pd.read_csv(path).drop(columns=['Unnamed: 0'])
    cluster_to_epi = {}
    for cluster in df.cluster.unique():
        cluster_to_epi[cluster] = df[df.cluster == cluster]
    return cluster_to_epi


if __name__ == "__main__":
    alpha_cluster = 29
    beta_cluster = 2

    clean_beta_cm = pd.read_csv('../data/significant_clone_matrix_fisher_fmba_TRB_top_500k_wo_leaks.csv').drop(
        columns=['Unnamed: 0'])
    covid_clones_beta = clean_beta_cm.columns[1:]
    res_beta = seqs2hamming(covid_clones_beta, viz_method='drl')

    clean_alpha_cm = pd.read_csv('../data/significant_clone_matrix_fisher_fmba_TRA_top_500k_wo_leaks.csv').drop(
        columns=['Unnamed: 0'])
    covid_clones_alpha = clean_alpha_cm.columns[1:]
    res_alpha = seqs2hamming(covid_clones_alpha, viz_method='drl')

    vdjdb = pd.read_csv('../data/vdjdb.txt', sep='\t')
    ab_vdjdb = pd.read_csv('../data/vdjdb_full.txt', sep='\t')

    alpha_epi = set(get_epitopes_for_cluster(vdjdb, res_alpha, alpha_cluster, dist=1)['antigen.epitope'])
    beta_epi = set(get_epitopes_for_cluster(vdjdb, res_beta, beta_cluster, dist=1)['antigen.epitope'])

    intersected_epi = alpha_epi.intersection(beta_epi)

