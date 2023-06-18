import numpy as np
import pandas as pd
from scipy.stats import logistic, linregress
from sklearn.preprocessing import StandardScaler

pd.options.mode.chained_assignment = None  # default='warn'
epsilon = 0.000001


def normalize_usage_matrix_by_rows(matrix: pd.DataFrame):
    """
    A function that normalizes gene usages within one run.
    This is made to translate all the runs' data into vector normalized by 1.
    The given matrix can have any annotation columns, they will be thrown away,
    but the genes should be given in the following format: IGHV.-(.*)?
    :param matrix: the gene usage matrix, which is not preprocessed yet;
    :return: the normalized usage matrix without annotation columns
    """
    usage_matrix = matrix.copy()
    for column in usage_matrix.columns:
        if not column.startswith('TR'):
            usage_matrix = usage_matrix.drop(columns=[column])
    usage_matrix["sum"] = usage_matrix.sum(axis=1)
    for column in usage_matrix.columns:
        usage_matrix[column] = usage_matrix[column] / usage_matrix['sum']
    usage_matrix = usage_matrix.drop(columns=['sum'])
    return usage_matrix


def standardize_usage_matrix(usage_matrix: pd.DataFrame, annotation_columns):
    datasets = list(usage_matrix.project.unique())
    matrices = []
    for name in datasets:
        usage_matrix_name = usage_matrix[usage_matrix.project == name]
        annotation = usage_matrix_name[annotation_columns]
        v_genes = usage_matrix_name[[x for x in usage_matrix_name.columns if x.startswith('TR')]]
        v_gene_names = v_genes.columns
        v_genes = pd.DataFrame(data=StandardScaler().fit_transform(v_genes), columns=v_gene_names)
        matrices.append(pd.concat([annotation.reset_index(drop=True), v_genes], axis=1))
    return pd.concat(matrices).reset_index(drop=True)


def standardize_usage_matrix_log(usage_matrix: pd.DataFrame, annotation_columns):
    datasets = list(usage_matrix.project.unique())
    print(datasets)
    matrices = []
    for name in datasets:
        usage_matrix_name = usage_matrix[usage_matrix['project'] == name]
        annotation = usage_matrix_name[annotation_columns]
        v_genes = usage_matrix_name[[x for x in usage_matrix_name.columns if x.startswith('TR')]]
        v_genes = v_genes.apply(lambda x: np.log(x + epsilon))
        v_gene_names = v_genes.columns
        v_genes = pd.DataFrame(data=StandardScaler().fit_transform(v_genes), columns=v_gene_names)
        matrices.append(pd.concat([annotation.reset_index(drop=True), v_genes], axis=1))
    cur_usage_matrix = pd.concat(matrices).reset_index(drop=True)
    return cur_usage_matrix


def get_sigmoid_parameters_for_v_gene(raw_freq: pd.Series, stand_freq: pd.Series, v_gene):
    raw_freq += epsilon
    log_raw_freq = raw_freq.apply(lambda x: np.log((1 - x) / x))
    mdl = linregress(stand_freq, log_raw_freq)
    print(f'for v gene {v_gene} R-squared is {mdl.rvalue ** 2:.6f}, slope - {mdl.slope}, intercept - {mdl.intercept}')
    return mdl.slope, mdl.intercept


def standardize_usage_matrix_log_exp(usage_matrix: pd.DataFrame, annotation_columns):
    norm_usage_matrix = normalize_usage_matrix_by_rows(usage_matrix)
    log_stand_usage_matrix = standardize_usage_matrix_log(usage_matrix, annotation_columns)
    for name in log_stand_usage_matrix.columns:
        if 'TR' in name:
            a = 1
            b = norm_usage_matrix[name].mean()
            log_stand_usage_matrix[name] = log_stand_usage_matrix[name].apply(lambda x: 2 * b * logistic.cdf(x))
    return log_stand_usage_matrix


def prepare_usage_matrix(usage_matrix: pd.DataFrame,
                         annotation_columns=['run', 'project', 'covid', 'hla', 'number_of_clonotypes'],
                         standardize_method=None):
    norm_usage_matrix = normalize_usage_matrix_by_rows(usage_matrix)
    normalized_usage_matrix = usage_matrix[annotation_columns]
    for col in norm_usage_matrix.columns:
        normalized_usage_matrix.loc[:, col] = norm_usage_matrix[col]
    if standardize_method is not None:
        normalized_usage_matrix = standardize_method(normalized_usage_matrix, annotation_columns)
    return normalized_usage_matrix


def create_usage_matrices_for_fmba_beta():
    suffix = 'fmba_TRB'
    usage_matrix = pd.read_csv(f'data/usage_matrix_{suffix}.csv').drop(columns=['Unnamed: 0'])
    norm_um = prepare_usage_matrix(usage_matrix)
    norm_um.to_csv(f'data/normalized_usage_matrix_{suffix}.csv')
    prepare_usage_matrix(usage_matrix, standardize_method=standardize_usage_matrix_log_exp).to_csv(
        f'data/standardized_usage_matrix_{suffix}_wo_test_runs.csv')


def create_usage_matrices_for_fmba_alpha():
    suffix = 'fmba_TRA'
    usage_matrix = pd.read_csv(f'data/usage_matrix_{suffix}.csv').drop(columns=['Unnamed: 0'])
    norm_um = prepare_usage_matrix(usage_matrix)
    norm_um.to_csv(f'data/normalized_usage_matrix_{suffix}.csv')
    prepare_usage_matrix(usage_matrix, standardize_method=standardize_usage_matrix_log_exp).to_csv(
        f'data/standardized_usage_matrix_{suffix}_wo_test_runs.csv')


def create_usage_matrices_for_adaptive():
    suffix = 'adaptive'
    usage_matrix = pd.read_csv(f'data/usage_matrix_{suffix}.csv').drop(columns=['Unnamed: 0']).dropna(axis=1)
    norm_um = prepare_usage_matrix(usage_matrix, annotation_columns=['run', 'project', 'covid'])
    norm_um.to_csv(f'data/normalized_usage_matrix_{suffix}.csv')
    prepare_usage_matrix(usage_matrix, standardize_method=standardize_usage_matrix_log_exp,
        annotation_columns=['run', 'project', 'covid']).to_csv(
        f'data/standardized_usage_matrix_{suffix}.csv')


def create_usage_matrices_for_joint():
    suffix = 'joint_new'
    usage_matrix = pd.read_csv(f'../data/usage_matrix_{suffix}.csv').drop(columns=['Unnamed: 0']).dropna(axis=1)
    norm_um = prepare_usage_matrix(usage_matrix, annotation_columns=['run', 'project', 'covid'])
    norm_um.to_csv(f'../data/normalized_usage_matrix_{suffix}.csv')
    prepare_usage_matrix(usage_matrix, standardize_method=standardize_usage_matrix_log_exp,
        annotation_columns=['run', 'project', 'covid']).to_csv(
        f'../data/standardized_log_exp_usage_matrix_{suffix}.csv')


def create_usage_matrices_for_functional():
    for suffix in ['functional', 'nonfunctional']:
        usage_matrix = pd.read_csv(f'../data/usage_matrix_{suffix}.csv').drop(columns=['Unnamed: 0'])
        norm_um = prepare_usage_matrix(usage_matrix)
        norm_um.to_csv(f'../data/normalized_usage_matrix_{suffix}.csv')
        prepare_usage_matrix(usage_matrix, standardize_method=standardize_usage_matrix_log_exp).to_csv(
            f'../data/standardized_log_exp_usage_matrix_{suffix}.csv')


def create_joint_TRB_adaptive_fmba_um():
    suffix = 'joint'
    adaptive_stand_um = pd.read_csv('data/usage_matrix_adaptive.csv').drop(columns=['Unnamed: 0'])
    fmba_stand_um = pd.read_csv('data/usage_matrix_fmba_TRB.csv').drop(columns=['Unnamed: 0'])
    joint_um = pd.concat([fmba_stand_um, adaptive_stand_um])
    joint_um = joint_um[['run', 'project', 'covid'] + [x for x in joint_um.columns if x.startswith('TRB')]].fillna(0)
    norm_um = prepare_usage_matrix(joint_um, annotation_columns=['run', 'project', 'covid'])
    norm_um.to_csv(f'data/normalized_usage_matrix_{suffix}.csv')
    prepare_usage_matrix(joint_um, standardize_method=standardize_usage_matrix_log_exp,
                         annotation_columns=['run', 'project', 'covid']).to_csv(
        f'data/standardized_usage_matrix_{suffix}.csv')


if __name__ == "__main__":
    if 'snakemake' in globals():
        if snakemake.params.gene == 'TRB':
            if snakemake.params.platform == 'fmba':
                create_usage_matrices_for_fmba_beta()
            elif snakemake.params.platform == 'adaptive':
                create_usage_matrices_for_adaptive()
            elif snakemake.params.platform == 'joint':
                create_joint_TRB_adaptive_fmba_um()
        if snakemake.params.gene == 'TRA':
            if snakemake.params.platform == 'fmba':
                create_usage_matrices_for_fmba_alpha()
