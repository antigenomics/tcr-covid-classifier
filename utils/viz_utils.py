import math
from math import sqrt, pi

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde as kde
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import Normalizer

from utils.clustering_utils import get_most_frequent_cluster_by_vdjdb_occurence
from utils.data_utils import prepare_clonotype_matrix, prepare_run_column

import matplotlib.gridspec as gridspec

def plot_usage_matrix_pca(usage_matrix: pd.DataFrame, method=PCA, n_components=2, target=None, plot_gradient=False,
                          figsize=(12, 5), ax=None):
    def plot_results(pca_results, target, plot_name):
        pca_results[plot_name] = target
        if plot_gradient:
            if len(target.unique()) == 1:
                print('Bad gene:(')
                return
            # ax = plt.axes()
            if n_components > 2:
                sns.pairplot(pca_results, plot_kws=dict(
                    hue=target,
                    palette=mpl.cm.viridis
                ))
            else:
                if ax is None:
                    plt.figure(figsize=figsize, dpi=80)
                    sc = plt.scatter(pca_results['PC1'], pca_results['PC2'], c=target, vmin=min(target),
                                     vmax=max(target))
                    plt.colorbar(sc)
                else:
                    sc = ax.scatter(pca_results['PC1'], pca_results['PC2'], c=target, vmin=min(target),
                                    vmax=max(target))
                    plt.colorbar(sc, ax=ax)
        else:
            if n_components > 2:
                sns.pairplot(pca_results, hue=plot_name)
            else:
                sns.scatterplot(x='PC1', y='PC2', data=pca_results, hue=plot_name, ax=ax)
                sns.move_legend(ax, "upper right", bbox_to_anchor=(1, 1))
        if ax is None:
            plt.show()
        ax.set_xlabel('Component 1')
        ax.set_ylabel('Component 2')

    pca = method(n_components=n_components, random_state=42)
    genes_names = [x for x in usage_matrix.columns if x.startswith('TR')]
    usage_matrix_genes = usage_matrix[genes_names]

    pca_results = pca.fit_transform(usage_matrix_genes)
    pca_results = pd.DataFrame(data=pca_results[:, :], columns=[f'PC{i + 1}' for i in range(n_components)])

    if target is None:
        for name in ['project', 'covid', 'hla']:
            plot_results(pca_results, usage_matrix[name], name)
    else:
        plot_results(pca_results, usage_matrix[target], target)
    if ax is None:
        plt.show()


def pca_results_appending(usage_matrix: pd.DataFrame):
    pca = PCA(n_components=2)
    genes_names = [x for x in usage_matrix.columns if x.startswith('TR')]
    usage_matrix_genes = usage_matrix[genes_names]

    pca_results = pca.fit_transform(usage_matrix_genes)
    pca_results = pd.DataFrame(data=pca_results[:, :], columns=['x', 'y'])
    usage_matrix['x'] = pca_results['x']
    usage_matrix['y'] = pca_results['y']
    return usage_matrix


def plot_multiple_pca(usage_matrix: pd.DataFrame, target):
    genes_names = [x for x in usage_matrix.columns if x.startswith('TR')]

    pca = PCA()
    components = pca.fit_transform(usage_matrix[genes_names])
    # labels = [f"PC {i + 1} ({var:.1f}%)"
    #     for i, var in enumerate(pca.explained_variance_ratio_ * 100)
    # ]
    labels = [i for i, var in enumerate(pca.explained_variance_ratio_ * 100)]

    df = pd.DataFrame({x: pd.Series(y) for x, y in zip(labels, components.T) if x < 6})
    df[target] = usage_matrix[target]
    sns.pairplot(df, hue=target)
    plt.show()


def plot_boxplots_for_usage_matrix(usage_matrix, by, filename=None, layout=(10, 8), figsize=(70, 50), log_scale=False,
                                   ylim=None):
    plt.clf()
    plt.rcParams["figure.figsize"] = figsize
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    col = []
    for column_name in usage_matrix.columns:
        if column_name.startswith('TR') or column_name == by:
            col.append(column_name)
    important_matrix = usage_matrix[col]
    if log_scale:
        for name in col:
            if name.startswith('TR'):
                important_matrix[name] = important_matrix[name].apply(lambda x: math.log(x + 1))
    # for column_name in important_matrix.columns:
    #     if column_name.startswith('IGH'):
    #         important_matrix[column_name] = important_matrix[column_name] / important_matrix[column_name].max()
    plot = important_matrix[col].boxplot(by=by, layout=layout, return_type='axes', ax=ax, rot=90)

    if ylim is not None:
        for ax_i in plot:
            ax_i.set_ylim(ylim)
    if filename is not None:
        fig.savefig(filename)


def _plot_feature_importance(importance, names, model_type, k_best=20, ax=None):
    # Create arrays from feature importance and feature names
    feature_importance = np.array(importance)
    feature_names = np.array(names)

    # Create a DataFrame using a Dictionary
    data = {'feature_names': feature_names, 'feature_importance': feature_importance}
    fi_df = pd.DataFrame(data)

    # Sort the DataFrame in order decreasing feature importance
    fi_df.sort_values(by=['feature_importance'], ascending=False, inplace=True)
    fi_df = fi_df.head(k_best)

    # Define size of bar plot
    if ax is None:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot()
    # Plot Searborn bar chart
    sns.barplot(x=fi_df['feature_names'], y=fi_df['feature_importance'], ax=ax, orient='v')
    # Add chart labels
    ax.set_title(model_type)
    ax.set_xlabel('FEATURE IMPORTANCE')
    ax.set_ylabel('FEATURE NAMES')
    ax.tick_params(labelrotation=90, axis='x')


def plot_feature_importances(model, feature_names, model_type, k_best, ax=None):
    feature_importances = model.feature_importances_
    _plot_feature_importance(feature_importances, feature_names, model_type, k_best, ax)


def plot_v_usage_hist(usage_matrix: pd.DataFrame, target=None):
    fig, axes = plt.subplots(3, 3, figsize=[20, 20])
    axes = axes.flatten()
    matrix = usage_matrix.groupby('project', as_index=False).mean()
    for i, proj in enumerate(matrix.project.unique()):
        cur_matrix = matrix[matrix.project == proj].T.iloc[1:, :]
        cur_matrix = cur_matrix.sort_values(by=cur_matrix.columns[0], ascending=False)
        cur_matrix = cur_matrix.reset_index().rename(columns={'index': 'v_gene', cur_matrix.columns[0]: 'usage'})
        _plot_feature_importance(cur_matrix['usage'], cur_matrix['v_gene'], f'v_gene for {proj} ', k_best=25,
                                 ax=axes[i])


def plot_feature_importance_rf(usage_matrix):
    usage_matrix_dropped = usage_matrix.drop(columns=['run', 'project', 'covid', 'hla'])
    X = usage_matrix_dropped
    X = Normalizer().fit_transform(X)
    y = usage_matrix['covid']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.15, random_state=42)
    rf = RandomForestClassifier(random_state=42)
    rf.fit(X_train, y_train)
    print(f'Accuracy score is {accuracy_score(rf.predict(X_test), y_test)}')
    _plot_feature_importance(rf.feature_importances_, usage_matrix_dropped.columns, 'RANDOM FOREST', k_best=25)


def plot_cluster_map(usage_matrix, color_by='platform', ax=None):
    v_gene_columns = [x for x in usage_matrix.columns if 'TR' in x]

    labels = usage_matrix[color_by]
    pal = sns.cubehelix_palette(labels.unique().size,
                                light=.9, dark=.1, reverse=True,
                                start=1, rot=-2)
    lut = dict(zip(map(str, labels.unique()), pal))
    label_cols = labels.map(lut)

    g = sns.clustermap(usage_matrix[v_gene_columns],
                       xticklabels=True,
                       cmap='vlag',
                       # figsize=(20, 25),
                       row_colors=label_cols)

    for label in labels.unique():
        g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                                label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc="center", ncol=6)
    return g


def plot_clustermap_axes_based(usage_matrix, genes=['TRBV28', 'TRBV4-3', 'TRBV6-2'], ax=None):
    v_gene_columns = [x for x in usage_matrix.columns if 'TR' in x]
    Z = linkage(usage_matrix[genes], 'complete')
    labels = fcluster(Z, 0.00001, criterion='distance')
    labels_order = np.argsort(labels)
    needed_data = usage_matrix[genes]
    for gene in genes:
        needed_data[gene] = (usage_matrix[gene] - usage_matrix[gene].min()) / (usage_matrix[gene].max() - usage_matrix[gene].min())
    sns.heatmap(needed_data.loc[labels_order, :][genes].T, ax=ax, cbar_kws={'label': 'V gene usage in a sample'},
                cmap='vlag')
    for k, gene in enumerate(genes):
        zscores = zscore(usage_matrix.loc[labels_order, :][gene])
        sns.lineplot(x=[i for i in range(len(usage_matrix[gene]))],
                     y=-(zscores - zscores.min()) / (zscores.max() - zscores.min()) + 1 + k, ax=ax, color='black', linewidth=0.5)
    ax.get_xaxis().set_visible(False)


from scipy.stats import mannwhitneyu, linregress
class SeabornFig2Grid():

    def __init__(self, seaborngrid, fig,  subplot_spec):
        self.fig = fig
        self.sg = seaborngrid
        self.subplot = subplot_spec
        if isinstance(self.sg, sns.axisgrid.FacetGrid) or \
            isinstance(self.sg, sns.axisgrid.PairGrid):
            self._movegrid()
        elif isinstance(self.sg, sns.axisgrid.JointGrid):
            self._movejointgrid()
        self._finalize()

    def _movegrid(self):
        """ Move PairGrid or Facetgrid """
        self._resize()
        n = self.sg.axes.shape[0]
        m = self.sg.axes.shape[1]
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n,m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i,j], self.subgrid[i,j])

    def _movejointgrid(self):
        """ Move Jointgrid """
        h= self.sg.ax_joint.get_position().height
        h2= self.sg.ax_marg_x.get_position().height
        r = int(np.round(h/h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(r+1,r+1, subplot_spec=self.subplot)

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax, gs):
        #https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure=self.fig
        self.fig.axes.append(ax)
        self.fig.add_axes(ax)
        ax._subplotspec = gs
        ax.set_position(gs.get_position(self.fig))
        ax.set_subplotspec(gs)

    def _finalize(self):
        plt.close(self.sg.fig)
        self.fig.canvas.mpl_connect("resize_event", self._resize)
        self.fig.canvas.draw()

    def _resize(self, evt=None):
        self.sg.fig.set_size_inches(self.fig.get_size_inches())


def significant_clones_distribution(significant_clonotype_matrix, run_to_number_of_clones, desc, data_description,
                                    by='covid', reg_ax=None, hist_ax=None, fit_reg=True):
    for df in [significant_clonotype_matrix, run_to_number_of_clones, desc]:
        prepare_run_column(df)
    sign_clone_count = significant_clonotype_matrix.shape[1] - 1
    significant_clonotype_matrix['sum'] = significant_clonotype_matrix.sum(axis=1)

    significant_clonotype_matrix = significant_clonotype_matrix.merge(desc[['run', by]]).merge(run_to_number_of_clones)
    healthy = significant_clonotype_matrix[
        significant_clonotype_matrix[by] == ('-' if by == 'cmv' else 'healthy')].drop(columns=[by])
    ill = significant_clonotype_matrix[significant_clonotype_matrix[by] != ('-' if by == 'cmv' else 'healthy')].drop(
        columns=[by])
    if reg_ax is None and hist_ax is None:
        f, (reg_ax, hist_ax) = plt.subplots(1, 2, figsize=(15, 6))

    if hist_ax is not None:
        print('updated version')
        bins_count = 10 if len(healthy) < 300 else 25

        healthy_percents = pd.DataFrame(data={'percent': healthy['sum'] / healthy['number_of_clones']})
        healthy_percents['status'] = 'healthy'

        ill_percents = pd.DataFrame(data={'percent': ill['sum'] / ill['number_of_clones']})
        ill_percents['status'] = by

        sns.histplot(data=pd.concat([healthy_percents, ill_percents]),
                     x='percent',
                     kde=True,
                     alpha=0.5,
                     stat="density",
                     ax=hist_ax,
                     hue='status',
                     bins=bins_count
                     )

        p = '{:.3g}'.format(mannwhitneyu(healthy_percents.percent, ill_percents.percent)[1])

        hist_ax.set_ylabel(f'fraction of samples')
        hist_ax.set_xlabel(f'fraction of {by} associated clones in a sample')
        hist_ax.set_title(f'{sign_clone_count} {data_description}, M-W P={p}')

    if reg_ax is not None:
        slope_h, intercept_h, r_h, p_h, sterr_h = linregress(x=healthy['number_of_clones'],
                                                             y=healthy['sum'])
        equation_healthy = 'healthy (y = ' + str(round(slope_h, 3)) + 'x{0:+g}'.format(round(intercept_h, 3)) + ')'
        print(slope_h, intercept_h, r_h, p_h, sterr_h)
        sns.regplot(x='number_of_clones',
                    y='sum',
                    data=healthy,
                    label=equation_healthy,
                    scatter_kws={'alpha': 0.5},
                    ax=reg_ax,
                    fit_reg=fit_reg)

        slope_c, intercept_c, r_c, p_c, sterr_c = linregress(x=ill['number_of_clones'],
                                                             y=ill['sum'])
        print(slope_c, intercept_c, r_c, p_c, sterr_c)
        equation_covid = 'covid (y = ' + str(round(slope_c, 3)) + 'x + ' + str(round(intercept_c, 3)) + ')'
        sns.regplot(x='number_of_clones',
                    y='sum',
                    data=ill,
                    label=equation_covid,
                    scatter_kws={'alpha': 0.5},
                    ax=reg_ax,
                    fit_reg=fit_reg)
        reg_ax.set_xlabel('count of unique TCRs')
        reg_ax.set_ylabel(f'count of {by} associated clones')
        reg_ax.set_title(f'{by} associated clones for {data_description}')
        reg_ax.legend()


def plot_results_for_hla_class(hla_class, hla_keys, percents, features_counts, best_scores, errors=None,
                               val_scores=None):
    keys_for_hla = [x for x in hla_keys if hla_class in x and x != 'A*69']
    keys_for_hla.sort(reverse=True, key=lambda x: percents[x])
    if val_scores is None:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 6))
    else:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(20, 6))
    if errors is None:
        ax1.barh(range(len(keys_for_hla)), [percents[x] for x in keys_for_hla], align='center')
    else:
        ax1.barh(range(len(keys_for_hla)), [percents[x] for x in keys_for_hla], align='center')

    ax1.set_yticks(range(len(keys_for_hla)), keys_for_hla)
    ax1.set_title(f'HLA-{hla_class} alleles frequencies in FMBA data')
    ax1.set_xlabel('% of data with allele')
    ax1.invert_yaxis()

    ax2.barh(range(len(keys_for_hla)), [features_counts[x] for x in keys_for_hla], align='center')
    ax2.set_yticks(range(len(keys_for_hla)), keys_for_hla)
    ax2.set_title(f'Count of clones associated with HLA-{hla_class} alleles')
    ax2.set_xlabel('# of HLA-associated clonotypes')
    ax2.invert_yaxis()

    ax3.barh(range(len(keys_for_hla)), [best_scores[x] for x in keys_for_hla], xerr=[errors[x] for x in keys_for_hla],
             align='center')
    ax3.set_yticks(range(len(keys_for_hla)), keys_for_hla)
    ax3.set_title(f'F1-scores for prediction of HLA-{hla_class} alleles')
    ax3.set_xlabel('f1-score')
    ax3.invert_yaxis()

    if val_scores is not None:
        ax4.barh(range(len(keys_for_hla)), [val_scores[x] for x in keys_for_hla], align='center')
        ax4.set_yticks(range(len(keys_for_hla)), keys_for_hla)
        ax4.set_title(f'F1-scores for prediction of HLA-{hla_class} alleles in HIP')
        ax4.set_xlabel('f1-score')
        ax4.invert_yaxis()

    plt.show()


def plot_occurances_in_samples_distribution(clonotype_matrix_path, desc_path, bins=50):
    cm = prepare_clonotype_matrix(clonotype_matrix_path, make_bool_features=True)
    cm = cm.merge(prepare_clonotype_matrix(desc_path)[['run']])
    results = cm.drop(columns=['run']).sum(axis=0)
    results.hist(bins=bins)
    plt.xlabel('# of samples with clone')
    plt.ylabel('frequency')
    return results


def __create_data_for_platform(pgen_path, cm_path, run_to_number_of_clonotypes, make_bool_features=True):
    pgen = pd.read_csv(pgen_path, header=None, names=['clone', 'pgen'], error_bad_lines=False)
    sum_value = run_to_number_of_clonotypes.number_of_clones.sum()
    cm = prepare_clonotype_matrix(cm_path, make_bool_features=make_bool_features).drop(columns=['run'])
    preal = pd.DataFrame(cm.sum(axis=0) / sum_value).reset_index().rename(
        columns={'cdr3aa': 'clone', 0: 'preal', 'index': 'clone'})
    full_data = pgen.merge(preal)
    return full_data


def plot_generated_to_real_distribution(pgen_paths, cm_paths, desc, run_to_number_of_clones_path, labels=None,
                                        make_bool_features=True):
    if not isinstance(pgen_paths, list):
        pgen_paths = [pgen_paths]
        cm_paths = [cm_paths]
        labels = ['label']
    full_data = {}
    run_to_number_of_clonotypes = pd.read_csv(run_to_number_of_clones_path).set_index('run')
    for pgen_path, cm_path, label in zip(pgen_paths, cm_paths, labels):
        full_data[label] = __create_data_for_platform(pgen_path, cm_path, run_to_number_of_clonotypes,
                                                      make_bool_features)
    for label, data in full_data.items():
        sns.regplot(x='pgen', y='preal', data=data, label=label, scatter_kws={'alpha': 0.5}, logistic=True)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('log OLGA gen proba for clone')
    plt.ylabel('log Observed gen proba for clone')
    plt.legend()
    plt.title(f'Distribution of observed clone probabilities to generated ones, {desc}')
    plt.show()
    return full_data


def plot_olga_cleanup_data(probas, observed_threshold=1.2e-6, gen_threshold=3e-8, ax=None):
    if ax is None:
        fig, (ax) = plt.subplots(1, 1)
    sns.regplot(x='pgen', y='preal',
                data=probas[(probas['preal'] < observed_threshold) | (probas['pgen'] < gen_threshold)],
                label='rarely observed', scatter_kws={'alpha': 0.5}, fit_reg=False, ax=ax)
    sns.regplot(x='pgen', y='preal',
                data=probas[(probas['preal'] >= observed_threshold) & (probas['pgen'] >= gen_threshold)],
                label='often observed', scatter_kws={'alpha': 0.5}, fit_reg=False, ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('log OLGA generation probability')
    ax.set_ylabel('log Observed generation probability')
    ax.legend()


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


def plot_cluster_with_epitope_coloring(clustering_res, cluster_index, color_by='cluster', ax=None):
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


def plot_cooccurence_heatmap_with_epitopes_labeling(plotting_df, annot_df, ax=None, epitopes_count_threshold=0):
    if ax is None:
        if plotting_df.shape[0] > plotting_df.shape[1]:
            fig, ax = plt.subplots(figsize=(2.5, 20))
        else:
            fig, ax = plt.subplots(figsize=(20, 2.5))
    selected_cols = plotting_df.columns[(plotting_df.sum() > 0.1) & (annot_df.sum() >= epitopes_count_threshold)]
    sns.heatmap(
        data=plotting_df.loc[(plotting_df.sum(axis=1) > 0.1) & (annot_df.sum(axis=1) >= epitopes_count_threshold), :][
            selected_cols],
        ax=ax,
        linewidths=0.5,
        linecolor='black',
        cmap=sns.color_palette("coolwarm", as_cmap=True),
        annot=annot_df.loc[(plotting_df.sum(axis=1) > 0.1) & (annot_df.sum(axis=1) >= epitopes_count_threshold), :][
            selected_cols].applymap(lambda x: '.' if x == 0 else str(x)), fmt=''
    )
    ax.set_xlabel('α cluster index')
    ax.set_ylabel('β cluster index')
    # ax.set_title('α vs β cluster co-occurence matrix, co-occurence threshold=85%')


def plot_cooccurence_heatmap_with_epitopes_labeling_bubble(plotting_df, annot_df, fig=None, ax=None,
                                                           epitopes_count_threshold=0, corr_threshold_row=0.1,
                                                           corr_threshold_col=0.1, legend_x=1.05):
    if ax is None:
        if plotting_df.shape[0] > plotting_df.shape[1]:
            fig, ax = plt.subplots(figsize=(2.5, 20))
        else:
            fig, ax = plt.subplots(figsize=(20, 2.5))

    selected_cols = plotting_df.loc[:, (plotting_df.max() > corr_threshold_col) & (annot_df.sum() >= epitopes_count_threshold)].columns
    print(len(selected_cols))
    data_to_plot = \
    plotting_df.loc[(plotting_df.max(axis=1) > corr_threshold_row) & (annot_df.sum(axis=1) >= epitopes_count_threshold),selected_cols]
    bubble_sizes = \
    annot_df.loc[(plotting_df.max(axis=1) > corr_threshold_row) & (annot_df.sum(axis=1) >= epitopes_count_threshold), selected_cols]
    ylabels = np.array([x for x in data_to_plot.index])
    xlabels = np.array([x for x in data_to_plot.columns])

    N = len(ylabels)
    M = len(xlabels)
    x, y = np.meshgrid(np.arange(M), np.arange(N))
    s = np.array(bubble_sizes) + 1.5
    c = np.array(data_to_plot)

    R = s / s.max() / 2
    circles = [plt.Circle((j, i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
    col = PatchCollection(circles, array=c.flatten(), cmap="coolwarm")
    plot = ax.add_collection(col)

    ax.set(xticks=np.arange(M), yticks=np.arange(N),
           xticklabels=xlabels, yticklabels=ylabels)
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)
    ax.grid(which='minor')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='1%', pad=1)
    cbar = fig.colorbar(plot, cax=cax, orientation='vertical')
    # cbar.set_label('fraction of α-β co-occured pairs')
    cax.yaxis.set_ticks_position('left')
    cax.yaxis.set_label_position('left')
    cax.set_ylabel('fraction of α-β co-occured pairs')

    ax.set_xlabel('α cluster index')
    ax.set_ylabel('β cluster index')
    ax.tick_params(labelrotation=30)

    # line1 = plt.Line2D([], [], color="white", marker='o', markersize=8, markerfacecolor="slategray")
    # line2 = plt.Line2D([], [], color="white", marker='o', markersize=11, markerfacecolor="slategray")
    # line3 = plt.Line2D([], [], color="white", marker='o', markersize=14, markerfacecolor="slategray")
    # line4 = plt.Line2D([], [], color="white", marker='o', markersize=18, markerfacecolor="slategray")
    # line5 = plt.Line2D([], [], color="white", marker='o', markersize=21, markerfacecolor="slategray")
    # ax.legend((line2, line4),
    #           ('0', '1'),
    #           numpoints=1,
    #           loc=1,
    #           title='# epitopes',
    #           labelspacing=1.5,
    #           bbox_to_anchor=(legend_x, 1.15),
    #           frameon=False)


def plot_waterfall_by_column(data, proba_column, label_column, ax=None):
    results = data[[label_column, proba_column]].sort_values(by=proba_column)
    results = results.reset_index(drop=True).reset_index()
    results['probability'] = (results[proba_column] - 0.5) * 2
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 4))
    p = sns.barplot(data=results, x='index', y='probability', hue=label_column, ax=ax)
    p.set(xticklabels=[])
    p.set(xlabel=None)


def get_closest_delta(used_centers, cur_center):
    min_diff = 1000
    for center in used_centers:
        min_diff = min(min_diff, sqrt((center[0] - cur_center[0]) ** 2 + (center[1] - cur_center[1]) ** 2))
    return min_diff


def adjust_asin(x, y):
    if x > 0:
        if y > 0:
            return 0
        if y < 0:
            return pi / 2
    if x < 0:
        if y < 0:
            return pi
        return 3 * pi / 2


def create_epitope_name(epitope, vdjdb):
    mhc = list(vdjdb[vdjdb['antigen.epitope'] == epitope]['mhc.a'])[0]
    return mhc.split(':')[0].split('-')[1].replace('*', '') + '_' + epitope[:3]


def plot_clonotype_clustering_with_epitope_labeling(res, cluster_to_epi, vdjdb,
                                                    color_by='cluster',
                                                    cluster_size_threshold=0, dist_to_center=500,
                                                    center_diff_threshold=50, gene='TRB', ax=None,
                                                    global_zero_based=True,
                                                    species_to_plot=None):
    if ax is None:
        fig, (ax) = plt.subplots()
    plot_clusters_of_clonotypes(res, color_by=color_by, ax=ax)
    used_centers = set()
    cluster_num_to_epi_abbreviation_mapping = {}
    for k, v in cluster_to_epi.items():
        if v is not None:
            my_epi = v[~v['antigen.species'].str.contains('apiens')].reset_index(drop=True)
            if len(my_epi) > 0:
                for epitope, species in zip(my_epi['antigen.epitope'], my_epi['antigen.species']):
                    if species_to_plot is not None and species not in species_to_plot:
                        continue
                    cluster_center_x = list(res[res.cluster == k].x_mean)[0]
                    cluster_center_y = list(res[res.cluster == k].y_mean)[0]
                    cluster_size = list(res[res.cluster == k].cluster_size)[0]
                    if cluster_size > cluster_size_threshold:
                        dist_from_zero = sqrt(cluster_center_x ** 2 + cluster_center_y ** 2)
                        sinus_value = cluster_center_x / dist_from_zero
                        cosinus_value = cluster_center_y / dist_from_zero

                        if global_zero_based:
                            new_dist_to_plot = dist_to_center
                            if cluster_center_x < 0:
                                new_dist_to_plot += 100
                            new_point_x, new_point_y = new_dist_to_plot * sinus_value, new_dist_to_plot * cosinus_value
                            center_diff = get_closest_delta(used_centers, (new_point_x, new_point_y))

                            if center_diff < center_diff_threshold:
                                new_point_x = new_point_x - 200 * (1 if new_point_x < 0 else -1)
                                new_point_y = new_point_y - 200 * (1 if new_point_y < 0 else -1)
                                print(epitope)
                            used_centers.add((new_point_x, new_point_y))
                            ax.annotate(create_epitope_name(epitope, vdjdb),
                                        xy=(cluster_center_x, cluster_center_y),
                                        xytext=(new_point_x, new_point_y),
                                        bbox=dict(boxstyle='round,pad=0.2', fc='green' if species == 'SARS-CoV-2' else 'red', alpha=0.1),
                                        arrowprops={
                                            'arrowstyle': '->',
                                        })
                        else:
                            ax.annotate(create_epitope_name(epitope, vdjdb),
                                        xy=(cluster_center_x, cluster_center_y),
                                        xytext=(cluster_center_x + 2 * dist_to_center, cluster_center_y - dist_to_center),
                                        arrowprops={
                                            'arrowstyle': '->',
                                        })
                        cluster_num_to_epi_abbreviation_mapping[k] = create_epitope_name(epitope, vdjdb)
    ax.set_xlim(-600, 600)
    ax.set_ylim(-600, 600)
    ax.axis('off')
    return cluster_num_to_epi_abbreviation_mapping


def volcano_plot(clonotype_matrix_path, desc_path, pvals, healthy_col='covid', healthy_label='healthy',
                 platform='adaptive'):
    print('vvvvv')
    cm = prepare_clonotype_matrix(clonotype_matrix_path, make_bool_features=True).merge(
        prepare_run_column(pd.read_csv(desc_path)[['run', healthy_col, 'platform']])
    )
    cm = cm[cm.platform == platform]

    healthy_data = cm[cm[healthy_col] == healthy_label].drop_duplicates().set_index('run').drop(
        columns=['covid', 'platform']).T.reset_index().rename(
        columns={'index': 'clone'}).set_index('clone')
    healthy_data['count_of_ways_h'] = healthy_data.sum(axis=1)
    healthy_data = healthy_data.reset_index()[['count_of_ways_h', 'clone']]

    ill_data = cm[cm[healthy_col] != healthy_label].drop_duplicates().set_index('run').drop(
        columns=['covid', 'platform']).T.reset_index().rename(
        columns={'index': 'clone'}).set_index('clone')
    ill_data['count_of_ways_i'] = ill_data.sum(axis=1)
    ill_data = ill_data.reset_index()[['count_of_ways_i', 'clone']]
    print(healthy_data.merge(ill_data))

    df = pvals.merge(healthy_data).merge(ill_data)
    print(df)
    df['fold_change'] = df['count_of_ways_i'] / df['count_of_ways_h']
    df['log_fold_change'] = df['fold_change'].apply(lambda x: np.log2(x))
    df['logp'] = np.log2(df.pval + 1e-10)
    sns.scatterplot(data=df, x='log_fold_change', y='logp')


def make_different_colors(all_points, points):
    densObj = kde(all_points, bw_method=0.5)

    def makeColours(vals):
        colors = np.zeros((len(vals), 3))
        norm = Normalize(vmin=vals.min(), vmax=vals.max())
        colors = [cm.ScalarMappable(norm=norm, cmap='viridis').to_rgba(val) for val in vals]
        return colors

    return makeColours(densObj.evaluate(points))


def plot_volcano(data, pval_threshold=None, fold_change_threshold=None, selected_clones=None, pval_column='pval',
                 fold_change_column='log_fold_change', ax=None, pval_plotting_threshold=None):
    if ax is None:
        fig, ax = plt.subplots()
    if pval_threshold is not None and fold_change_threshold is not None:
        data['selected clone'] = data.apply(
            lambda row: True if row[fold_change_column] > fold_change_threshold and row[pval_column] < pval_threshold else False, axis=1)
    if selected_clones is not None:
        data['selected clone'] = data.clone.apply(lambda x: x in selected_clones)
    if pval_plotting_threshold is not None:
        data_to_plot = data[data[pval_column] < pval_plotting_threshold]
    else:
        data_to_plot = data
    plotting_data_background = data_to_plot[~data_to_plot['selected clone']][[fold_change_column, pval_column]].to_numpy().T
    plotting_data = data_to_plot[(data_to_plot['selected clone']) & (data_to_plot[fold_change_column] != np.inf)][
        [fold_change_column, pval_column]].to_numpy().T
    plotting_data[:, 1] = -np.log10(plotting_data[:, 1])

    plt.scatter(plotting_data_background[0], plotting_data_background[1], alpha=0.4, s=5,
                color=['grey' for _ in range(plotting_data_background.shape[1])])
    plt.scatter(plotting_data[0], plotting_data[1], alpha=0.4, s=5,
                color=make_different_colors(plotting_data, plotting_data))

    if fold_change_threshold is not None:
        ax.axvline(x=fold_change_threshold, linestyle='dashed', color='grey')
    if pval_threshold is not None:
        ax.axhline(y=pval_threshold, linestyle='dashed', color='grey')
    ax.set_xlabel('log fold change')
    ax.set_ylabel('p-value')

    ax.set_yscale('log')
    ax.invert_yaxis()


def plot_pandas_df_into_png(df, output_path=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
        # hide axes
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')
        fig.tight_layout()
    ax.table(cellText=df.values, colLabels=df.columns, loc='center')
    if output_path is not None:
        plt.savefig(output_path)
