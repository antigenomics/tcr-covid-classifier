import numpy as np
import pandas as pd
from tqdm import tqdm


def process_all_files(usage_matrix_path, save_path, method='top', count_of_clonotypes_to_choose=10000):
    um = pd.read_csv(usage_matrix_path).drop(columns=['Unnamed: 0', 'covid']).fillna(0)
    um = um[(um.run.str.contains('HIP')) | (um.run.str.contains('Keck'))]
    datasets_to_concat = []

    for run in tqdm(um['run'], total=len(um), desc='Extracting top clonotypes'):
        cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/hip_full/{run}', sep='\t')
        # cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/downsampled_new/{run}')
        cur_data = cur_data[['cdr3aa']].drop_duplicates()
        cur_data['count'] = 1
        datasets_to_concat.append(cur_data)
    full_data = pd.concat(datasets_to_concat)
    full_data = full_data[full_data.cdr3aa.str.isalpha()]

    if method == 'top':
        top = full_data.groupby(['cdr3aa'], as_index=False).count()
        top = top.sort_values(by=['count'], ascending=False).head(count_of_clonotypes_to_choose)

    elif method == 'random-roulette':
        top = full_data.groupby(['cdr3aa'], as_index=False).count()
        top = top.sample(n=count_of_clonotypes_to_choose, random_state=42, weights='count')

    elif method == 'random-uniform':
        top = full_data.groupby(['cdr3aa'], as_index=False).count()
        top = top.sample(n=count_of_clonotypes_to_choose, random_state=42)

    top.to_csv(save_path, index=False)


if __name__ == "__main__":
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=f'../data/most_used_clonotypes_hip_10000_wo_nt_uniform.csv',
                      method='random-uniform')
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=f'../data/most_used_clonotypes_hip_10000_wo_nt_roulette.csv',
                      method='random-roulette')
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=f'../data/most_used_clonotypes_hip_10000_raw.csv',
                      method='top')
