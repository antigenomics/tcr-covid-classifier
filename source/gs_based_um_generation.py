import math
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

from source.usage_matrix_preprocessing import prepare_usage_matrix

warnings.filterwarnings("ignore")


def get_golden_standard():
    olga = pd.read_csv('../data/olga/olga_generated_1e7_100.csv', sep='\t', header=None,
                       names=['cdr3nn', 'cdr3aa', 'v', 'j'])
    olga = olga.groupby(['v', 'j'], as_index=False).count().drop(columns=['cdr3aa']).rename(
        columns={'cdr3nn': 'gs_freq'})
    olga['gs_freq'] = olga['gs_freq'] / 1000000
    return olga


def evaluate_glm_model_poisson(data):
    formula = 'count ~ gs_freq + sample_size'
    model = smf.glm(formula=formula, data=data[['gs_freq', 'sample_size', 'count']],
                    family=sm.families.Poisson())
    fitted_model = model.fit()
    return fitted_model.predict(data)


def evaluate_glm_model_poisson_batch(data):
    formula = 'count ~ gs_freq + sample_size + C(v) + C(j)'
    model = smf.glm(formula=formula, data=data[['gs_freq', 'sample_size', 'count', 'v', 'j']],
                    family=sm.families.Poisson())
    fitted_model = model.fit()
    return fitted_model.predict(data)


def evaluate_glm_model_poisson_batch_with_project(data):
    formula = 'count ~ gs_freq + sample_size + C(v) + C(j) + C(project) + I(gs_freq * sample_size)'
    model = smf.glm(formula=formula, data=data[['gs_freq', 'sample_size', 'count', 'v', 'j', 'project']],
                    family=sm.families.Poisson())
    fitted_model = model.fit()
    return fitted_model.predict(data)


def evaluate_glm_model_binom(data):
    data['negative_count'] = data['sample_size'] - data['count']
    model = sm.GLM(data[['count', 'negative_count']], data[['gs_freq', 'sample_size']],
                   family=sm.families.Binomial())
    fitted_model = model.fit()
    print(fitted_model.summary())
    return fitted_model.predict(data[['gs_freq', 'sample_size']])


def evaluate_glm_model_binom_batch(data):
    data['negative_count'] = data['sample_size'] - data['count']
    model = sm.GLM(data[['count', 'negative_count']], data[['gs_freq', 'sample_size', 'count', 'j', 'v']],
                   family=sm.families.Binomial())
    fitted_model = model.fit()
    print(fitted_model.summary())
    return fitted_model.predict(data[['gs_freq', 'sample_size']])


def evaluate_glm_model_binom_batch_with_project(data):
    formula = 'count ~ gs_freq + sample_size + C(v) + C(j) + C(project) + I(gs_freq * sample_size)'
    model = smf.glm(formula=formula, data=data[['gs_freq', 'sample_size', 'count', 'v', 'j', 'project']],
                    family=sm.families.Binomial())
    fitted_model = model.fit()
    return fitted_model.predict(data)


def create_sample_size(vj_usage_matrix):
    sample_size_dict = defaultdict(int)
    for run, count in zip(vj_usage_matrix['run'], vj_usage_matrix['count']):
        sample_size_dict[run] += count
    vj_usage_matrix['sample_size'] = vj_usage_matrix['run'].apply(lambda x: sample_size_dict[x])
    return vj_usage_matrix, sample_size_dict


def adjust_vj_matrix(cur_vj_matrix, sample_size_dict, run_to_project_mapping, v, j, project):
    # 'v', 'j', 'count', 'run', 'project', 'sample_size', 'gs_freq'
    zero_runs = []
    zero_sample_sizes = []
    runs_with_vj_pair = set(cur_vj_matrix['run'])
    for run, sample_size in sample_size_dict.items():
        if run not in runs_with_vj_pair and run_to_project_mapping[run] == project:
            zero_runs.append(run)
            zero_sample_sizes.append(sample_size)
    zero_count = len(zero_runs)
    result = cur_vj_matrix.append(pd.DataFrame({
        'v': [v] * zero_count,
        'j': [j] * zero_count,
        'count': [0] * zero_count,
        'run': zero_runs,
        'project': [project] * zero_count,
        'sample_size': zero_sample_sizes
    })).reset_index().drop(columns=['index'])
    return result


def generate_v_gene_usage_vector(run, data):
    v_gene_usage_df = data[data['run'] == run]
    v_gene_to_usage_dict = {}
    for k, v in zip(v_gene_usage_df['v'], v_gene_usage_df['predicted_count']):
        v_gene_to_usage_dict[k] = v
    return v_gene_to_usage_dict


def fit_by_vj_pair(vj_usage_matrix, datasets, sample_size_dict, run_to_project_mapping, gs_data, evaluate_glm_model):
    datasets_to_concat = []
    for name in tqdm(datasets):
        cur_data = vj_usage_matrix[vj_usage_matrix['project'] == name]
        for v_gene in cur_data['v'].unique():
            v_gene_cur_data = cur_data[cur_data['v'] == v_gene]
            for j_gene in v_gene_cur_data['j'].unique():
                vj_cur_data = v_gene_cur_data[v_gene_cur_data['j'] == j_gene]
                vj_cur_data = adjust_vj_matrix(vj_cur_data, sample_size_dict, run_to_project_mapping, v_gene, j_gene,
                                               name)
                gs_vj_data = gs_data[(gs_data['v'] == v_gene) & (gs_data['j'] == j_gene)]

                if len(gs_vj_data) == 0:
                    vj_cur_data['predicted_count'] = 0
                    vj_cur_data['gs_freq'] = 0
                else:
                    vj_cur_data['gs_freq'] = list(gs_vj_data['gs_freq'])[0]
                    vj_cur_data['predicted_count'] = evaluate_glm_model(vj_cur_data)
                datasets_to_concat.append(vj_cur_data)
    full_data = pd.concat(datasets_to_concat)
    full_data = full_data.groupby(['v', 'run'], as_index=False).sum()
    return full_data


def fit_by_batch(vj_usage_matrix, datasets, sample_size_dict, run_to_project_mapping, gs_data, evaluate_glm_model):
    datasets_to_concat = []
    for name in tqdm(datasets):
        cur_data = vj_usage_matrix[vj_usage_matrix['project'] == name]
        merged_data = cur_data.merge(gs_data, on=['v', 'j'], how='left').fillna(0)
        merged_data['predicted_count'] = evaluate_glm_model(merged_data)
        datasets_to_concat.append(merged_data)
    full_data = pd.concat(datasets_to_concat)
    full_data = full_data.groupby(['v', 'run'], as_index=False).sum()
    return full_data


def fit_all_together(vj_usage_matrix, datasets, sample_size_dict, run_to_project_mapping, gs_data, evaluate_glm_model):
    merged_data = vj_usage_matrix.merge(gs_data, on=['v', 'j'], how='left').fillna(0)
    merged_data['predicted_count'] = evaluate_glm_model(merged_data)
    full_data = merged_data.groupby(['v', 'run'], as_index=False).sum()
    return full_data


def standardize_golden_standard_based(vj_usage_matrix: pd.DataFrame, usage_matrix: pd.DataFrame, evaluate_glm_model,
                                      grouping_way='vj',
                                      annotation_columns=['run', 'project', 'covid', 'hla', 'number_of_clonotypes']):
    datasets = list(vj_usage_matrix.project.unique())
    gs_data = get_golden_standard()
    vj_usage_matrix, sample_size_dict = create_sample_size(vj_usage_matrix)
    run_to_project_mapping = {}
    for run, proj in zip(usage_matrix['run'], usage_matrix['project']):
        run_to_project_mapping[run] = proj
    # нужно чтобы всегда были все возможные запуски для предсказания!
    if grouping_way == 'vj':
        full_data = fit_by_vj_pair(vj_usage_matrix, datasets, sample_size_dict, run_to_project_mapping, gs_data,
                                   evaluate_glm_model)
    elif grouping_way == 'batch':
        full_data = fit_by_batch(vj_usage_matrix, datasets, sample_size_dict, run_to_project_mapping, gs_data,
                                 evaluate_glm_model)
    elif grouping_way == 'none':
        full_data = fit_all_together(vj_usage_matrix, datasets, sample_size_dict, run_to_project_mapping, gs_data,
                                 evaluate_glm_model)

    all_v_genes = set()
    repertoires_list = []
    runs_series = full_data['run'].unique()
    for run in tqdm(runs_series):
        run_repertoire = generate_v_gene_usage_vector(run, full_data)
        all_v_genes.update(list(run_repertoire.keys()))
        repertoires_list.append(run_repertoire)

    full_usage_dict = {}
    for v_gene in tqdm(all_v_genes, desc='Making up usage dataframe'):
        v_gene_vector = []
        for repertoire_dict in repertoires_list:
            if v_gene in repertoire_dict:
                v_gene_vector.append(repertoire_dict[v_gene])
            else:
                v_gene_vector.append(0)
        full_usage_dict[v_gene] = v_gene_vector
    full_usage_dict['run'] = runs_series

    usage_df = pd.DataFrame.from_dict(full_usage_dict)

    return usage_df


if __name__ == "__main__":
    usage_matrix = pd.read_csv('../data/usage_matrix.csv').drop(columns=['Unnamed: 0'])
    vj_usage_matrix = pd.read_csv('../data/vj_usage_matrix.csv').drop(columns=['Unnamed: 0'])
    standardize_golden_standard_based(vj_usage_matrix, usage_matrix, evaluate_glm_model_poisson_batch_with_project,
                                      grouping_way='none').to_csv('../data/gs_usage_matrix_poisson_none_with_project.csv')
    standardize_golden_standard_based(vj_usage_matrix, usage_matrix, evaluate_glm_model_poisson_batch,
                                      grouping_way='none').to_csv('../data/gs_usage_matrix_poisson_none.csv')
    standardize_golden_standard_based(vj_usage_matrix, usage_matrix, evaluate_glm_model_poisson_batch,
                                      grouping_way='batch').to_csv('../data/gs_usage_matrix_poisson_batch.csv')
    # standardize_golden_standard_based(vj_usage_matrix, usage_matrix, evaluate_glm_model_binom).to_csv(
    #     '../data/gs_usage_matrix_binom.csv')
    standardize_golden_standard_based(vj_usage_matrix, usage_matrix, evaluate_glm_model_poisson,
                                      grouping_way='vj').to_csv('../data/gs_usage_matrix_poisson.csv')
