import datetime
import random
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from math import floor
from tqdm import tqdm

warnings.filterwarnings("ignore")

error_list = []


def update_counts_for_run(path_to_run, run_frequencies,
                          path_to_save=f'/projects/fmba_covid/1_data_links/downsampled_alpha'):
    run_name = path_to_run.split('/')[-1]
    print(f'started {run_name} at {datetime.datetime.now()}')
    try:
        raw_data = pd.read_csv(path_to_run, sep='\t')

        v_gene_names = set(run_frequencies.keys()).intersection(set(raw_data['v'].unique()))

        old_run_frequencies = run_frequencies.copy()
        run_frequencies = {x: y for (x, y) in old_run_frequencies.items() if x in v_gene_names}
        sum_freq = sum(run_frequencies.values())
        run_frequencies = {x: y / sum_freq for (x, y) in run_frequencies.items()}

        raw_data = raw_data[raw_data['v'].apply(lambda x: x in v_gene_names)]
        raw_data['old_freq'] = raw_data['freq'] / raw_data['freq'].sum()
        raw_data['old_count'] = raw_data['count']
        old_v_frequencies = {}
        for v in run_frequencies:
            old_v_frequencies[v] = raw_data[raw_data['v'] == v]['old_freq'].sum()

        full_count = raw_data['count'].sum()
        random_numbers = np.random.uniform(0, 1, full_count).tolist()
        random_numbers.sort()
        assert min(random_numbers) == random_numbers[0]
        assert max(random_numbers) == random_numbers[-1]

        raw_data['expected_freq'] = raw_data['v'].map(run_frequencies) * raw_data['old_freq'] / raw_data['v'].map(
            old_v_frequencies)
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
        error_list.append(path_to_run)

    print(f'ended {path_to_run} at {datetime.datetime.now()}')


def process_one_file(run):
    run_frequencies = transposed_um[[run]].to_dict()[run]
    if 'TCR' in run:
        dataset = 'adaptive_new'
    elif 'TR' in run:
        dataset = 'fmba_2021'
    else:
        dataset = 'hip_full'
    freqs_to_use = {x: y for (x, y) in run_frequencies.items() if 'TR' in x}
    update_counts_for_run(f'/projects/fmba_covid/1_data_links/{dataset}/{run}', freqs_to_use)


def process_all_files():
    runs = transposed_um.columns
    with Pool(72) as p:
        p.map(process_one_file, runs)


if __name__ == "__main__":
    print('started at', datetime.datetime.now())
    usage_matrix_path = '../data/standardized_log_exp_usage_matrix_by_v_gene_alpha_fmba.csv'
    um = pd.read_csv(usage_matrix_path).drop(columns=['Unnamed: 0']).fillna(0)
    um["sum"] = um.sum(axis=1, numeric_only=True)
    for column in um.columns:
        if column.startswith('TR'):
            um[column] = um[column] / um['sum']
    um = um.drop(columns=['sum'])
    transposed_um = um.set_index('run').transpose()
    process_all_files()
    print('finished at', datetime.datetime.now())
    print("errors in", error_list)