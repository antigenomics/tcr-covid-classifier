import datetime
import random
import warnings
from multiprocessing import Pool

import numpy as np
np.random.seed(42)
import pandas as pd
from math import floor
from tqdm import tqdm

warnings.filterwarnings("ignore")

error_list = []


def update_counts_for_run(path_to_run, run_frequencies_v, run_frequencies_j, path_to_save, count_of_clones_in_sample=None):
    run_name = path_to_run.split('/')[-1]
    # print(f'started {run_name} at {datetime.datetime.now()}')
    try:
        if '.clonotypes.' not in path_to_run:
            raw_data = pd.read_csv(path_to_run, sep='\t')
        else:
            raw_data = pd.read_csv(path_to_run)
        v_gene_names = set(run_frequencies_v.keys()).intersection(set(raw_data['v'].unique()))
        j_gene_names = set(run_frequencies_j.keys()).intersection(set(raw_data['j'].unique()))

        old_run_frequencies_v = run_frequencies_v.copy()
        old_run_frequencies_j = run_frequencies_j.copy()
        run_frequencies_v = {x: y for (x, y) in old_run_frequencies_v.items() if x in v_gene_names}
        run_frequencies_j = {x: y for (x, y) in old_run_frequencies_j.items() if x in j_gene_names}
        sum_freq_v = sum(run_frequencies_v.values())
        sum_freq_j = sum(run_frequencies_j.values())
        run_frequencies_v = {x: y / sum_freq_v for (x, y) in run_frequencies_v.items()}
        run_frequencies_j = {x: y / sum_freq_j for (x, y) in run_frequencies_j.items()}

        raw_data = raw_data[raw_data['v'].apply(lambda x: x in v_gene_names)]
        raw_data = raw_data[raw_data['j'].apply(lambda x: x in j_gene_names)]
        raw_data['old_freq'] = raw_data['freq'] / raw_data['freq'].sum()
        raw_data['old_count'] = raw_data['count']

        old_v_frequencies = {}
        for v in run_frequencies_v:
            old_v_frequencies[v] = raw_data[raw_data['v'] == v]['old_freq'].sum()
        old_j_frequencies = {}
        for j in run_frequencies_j:
            old_j_frequencies[j] = raw_data[raw_data['j'] == j]['old_freq'].sum()
        if count_of_clones_in_sample is None:
            full_count = raw_data['count'].sum()
        else:
            full_count = count_of_clones_in_sample
        random_numbers = np.random.uniform(0, 1, full_count).tolist()
        random_numbers.sort()
        assert min(random_numbers) == random_numbers[0]
        assert max(random_numbers) == random_numbers[-1]

        raw_data['expected_freq'] = raw_data['old_freq'] * (raw_data['v'].map(run_frequencies_v) / raw_data['v'].map(
            old_v_frequencies)) * (raw_data['j'].map(run_frequencies_j) / raw_data['j'].map(old_j_frequencies))
        raw_data['expected_freq'] = raw_data['expected_freq'] / raw_data['expected_freq'].sum()
        assert abs(raw_data['expected_freq'].sum() - 1) < 0.0005
        raw_data = raw_data.sort_values(by=['expected_freq']).reset_index(drop=True)
        sorted_frequencies = raw_data['expected_freq'].tolist()
        assert min(sorted_frequencies) == sorted_frequencies[0]
        assert max(sorted_frequencies) == sorted_frequencies[-1]

        for i in range(1, len(sorted_frequencies)):
            sorted_frequencies[i] += sorted_frequencies[i - 1]
        assert abs(sorted_frequencies[-1] - 1) < 0.0005
        raw_data['prefix_sum'] = pd.Series(sorted_frequencies)

        new_counts = [0 for _ in range(len(sorted_frequencies))]

        random_pointer = 0
        freq_pointer = 0
        while random_pointer < full_count:
            if random_numbers[random_pointer] <= sorted_frequencies[freq_pointer]:
                new_counts[freq_pointer] += 1
                random_pointer += 1
            else:
                freq_pointer += 1

        assert sum(new_counts) == full_count

        raw_data['count'] = pd.Series(new_counts)
        raw_data['freq'] = raw_data['count'] / raw_data['count'].sum()
        raw_data[raw_data['count'] > 0].loc[::-1].reset_index(drop=True).to_csv(path_to_save + f'/{run_name}',
                                                                                    index=False)
    except Exception as e:
        print(e, '!!')
        error_list.append(path_to_run)



def fmba_process_run(run, freqs_to_use_v, freqs_to_use_j, save_path):
    update_counts_for_run(f'{raw_data_path}/fmba/{run}',
                          run_frequencies_v=freqs_to_use_v,
                          run_frequencies_j=freqs_to_use_j,
                          path_to_save=save_path)


def adaptive_process_run(run, freqs_to_use_v, freqs_to_use_j, save_path, count_of_clones_in_sample=None):
    if 'HIP' in run or 'Keck' in run:
        update_counts_for_run(f'{raw_data_path}/hip_full/{run}',
                              run_frequencies_v=freqs_to_use_v,
                              run_frequencies_j=freqs_to_use_j,
                              path_to_save=save_path,
                              count_of_clones_in_sample=count_of_clones_in_sample)
    else:
        update_counts_for_run(f'{raw_data_path}/adaptive_new/{run}',
                              run_frequencies_v=freqs_to_use_v,
                              run_frequencies_j=freqs_to_use_j,
                              path_to_save=save_path,
                              count_of_clones_in_sample=count_of_clones_in_sample)


def process_one_file(run):
    run_frequencies_v = transposed_um_v[[run]].to_dict()[run]
    freqs_to_use_v = {x: y for (x, y) in run_frequencies_v.items() if 'TR' in x}

    run_frequencies_j = transposed_um_j[[run]].to_dict()[run]
    freqs_to_use_j = {x: y for (x, y) in run_frequencies_j.items() if 'TR' in x}
    if platform == 'fmba':
        fmba_process_run(run, freqs_to_use_v, freqs_to_use_j, f'{raw_data_path}/downsampled_{platform}_{gene}')
    elif platform == 'adaptive':
        adaptive_process_run(run=run,
                             freqs_to_use_v=freqs_to_use_v,
                             freqs_to_use_j=freqs_to_use_j,
                             save_path=f'{raw_data_path}/downsampled_{platform}_{gene}',
                             count_of_clones_in_sample=50000)
    elif platform == 'joint':
        if '.clonotypes.TRB' in run:
            fmba_process_run(run, freqs_to_use_v, freqs_to_use_j, f'{raw_data_path}/downsampled_joint')
        else:
            adaptive_process_run(run, freqs_to_use_v, freqs_to_use_j, f'{raw_data_path}/downsampled_joint')
    elif platform == 'joint_50k':
        if '.clonotypes.TRB' in run:
            fmba_process_run(run, freqs_to_use_v, freqs_to_use_j, f'{raw_data_path}/downsampled_joint_50k')
        else:
            adaptive_process_run(run, freqs_to_use_v, freqs_to_use_j, f'{raw_data_path}/downsampled_joint_50k',
                                 count_of_clones_in_sample=50000)


def process_all_files():
    runs = transposed_um_v.columns
    with Pool(snakemake.threads) as p:
        p.map(process_one_file, runs)


if __name__ == "__main__":
    print('started at', datetime.datetime.now())

    if 'snakemake' in globals():
        usage_matrix_path_v = snakemake.input[0]
        usage_matrix_path_j = snakemake.input[1]
        platform = snakemake.params.platform
        gene = snakemake.params.gene
        raw_data_path = snakemake.config['all_raw_data_path']
        save_path = snakemake.output[0]

        import os
        os.mkdir(save_path)

        um_v = pd.read_csv(usage_matrix_path_v).drop(columns=['Unnamed: 0']).fillna(0)
        um_v["sum"] = um_v.sum(axis=1, numeric_only=True)
        for column in um_v.columns:
            if column.startswith('TR'):
                um_v[column] = um_v[column] / um_v['sum']
        um_v = um_v.drop(columns=['sum'])

        um_j = pd.read_csv(usage_matrix_path_j).drop(columns=['Unnamed: 0']).fillna(0)
        print(um_j)
        um_j["sum"] = um_j.sum(axis=1, numeric_only=True)
        for column in um_j.columns:
            if column.startswith('TR'):
                um_j[column] = um_j[column] / um_j['sum']
        um_j = um_j.drop(columns=['sum'])

        transposed_um_v = um_v.set_index('run').transpose()
        transposed_um_j = um_j.set_index('run').transpose()
        process_all_files()

    print('finished at', datetime.datetime.now())
    print("errors in", error_list)
