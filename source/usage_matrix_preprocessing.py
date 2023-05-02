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
        # norm_usage_matrix = normalize_usage_matrix_by_rows(cur_usage_matrix)
        # normalized_usage_matrix = cur_usage_matrix[annotation_columns]
        # for col in norm_usage_matrix.columns:
        #     normalized_usage_matrix.loc[:, col] = norm_usage_matrix[col]
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
            # a, b = get_sigmoid_parameters_for_v_gene(raw_freq=norm_usage_matrix[name],
            #                                          stand_freq=log_stand_usage_matrix[name],
            #                                          v_gene=name)
            a = 1
            b = norm_usage_matrix[name].mean()
            print(f'b for {name} is {b}')
            log_stand_usage_matrix[name] = log_stand_usage_matrix[name].apply(lambda x: 2 * b * logistic.cdf(x))
    return log_stand_usage_matrix


def prepare_usage_matrix(usage_matrix: pd.DataFrame,
                         annotation_columns=['run', 'project', 'covid', 'hla', 'number_of_clonotypes'],
                         standardize_method=None):
    print(usage_matrix)
    norm_usage_matrix = normalize_usage_matrix_by_rows(usage_matrix)
    print(norm_usage_matrix)
    normalized_usage_matrix = usage_matrix[annotation_columns]
    for col in norm_usage_matrix.columns:
        normalized_usage_matrix.loc[:, col] = norm_usage_matrix[col]
    if standardize_method is not None:
        normalized_usage_matrix = standardize_method(normalized_usage_matrix, annotation_columns)
    return normalized_usage_matrix


def create_usage_matrices_for_fmba():
    suffix = ''
    usage_matrix = pd.read_csv('../data/usage_matrix' + suffix + '.csv').drop(columns=['Unnamed: 0'])
    norm_um = prepare_usage_matrix(usage_matrix)
    norm_um.to_csv('../data/normalized_usage_matrix.csv')
    prepare_usage_matrix(usage_matrix, standardize_method=standardize_usage_matrix_log_exp).to_csv(
        '../data/standardized_log_exp_usage_matrix_by_' + suffix + '.csv')


def create_usage_matrices_for_fmba_alpha():
    suffix = 'alpha_fmba'
    usage_matrix = pd.read_csv(f'../data/usage_matrix_{suffix}.csv').drop(columns=['Unnamed: 0'])
    norm_um = prepare_usage_matrix(usage_matrix)
    norm_um.to_csv(f'../data/normalized_usage_matrix_{suffix}.csv')
    prepare_usage_matrix(usage_matrix, standardize_method=standardize_usage_matrix_log_exp).to_csv(
        f'../data/standardized_log_exp_usage_matrix_by_{suffix}.csv')


def create_usage_matrices_for_adaptive():
    suffix = 'adaptive'
    usage_matrix = pd.read_csv(f'../data/usage_matrix_{suffix}.csv').drop(columns=['Unnamed: 0']).dropna(axis=1)
    norm_um = prepare_usage_matrix(usage_matrix, annotation_columns=['run', 'project', 'covid'])
    norm_um.to_csv(f'../data/normalized_usage_matrix_{suffix}.csv')
    prepare_usage_matrix(usage_matrix, standardize_method=standardize_usage_matrix_log_exp,
        annotation_columns=['run', 'project', 'covid']).to_csv(
        f'../data/standardized_log_exp_usage_matrix_by_{suffix}.csv')


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


if __name__ == "__main__":
    # create_usage_matrices_for_fmba_alpha()
    create_usage_matrices_for_joint()
    # create_usage_matrices_for_adaptive()
    # create_usage_matrices_for_functional()
