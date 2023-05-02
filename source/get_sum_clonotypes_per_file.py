from multiprocessing import Manager, Pool

import numpy as np
import pandas as pd
from tqdm import tqdm
run_to_clone_number = Manager().dict()


def process_one_file(run):
    try:
        cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/downsampled_alpha/{run}')
        # cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/hip_full/{run}', sep='\t')
        run_to_clone_number[run] = cur_data.shape[0]
        print(len(run_to_clone_number))
    except Exception as e:
        pass


def process_all_files():
    runs = list(um['run'])
    print(len(runs))
    with Pool(80) as p:
        p.map(process_one_file, runs)
    data_dict = {'run': [], 'number_of_clones': []}
    for x, y in run_to_clone_number.items():
        data_dict['run'].append(x)
        data_dict['number_of_clones'].append(y)
    pd.DataFrame(data=data_dict).to_csv(f'../data/run_to_number_of_clonotypes_alpha.csv', index=False)


if __name__ == "__main__":
    um = pd.read_csv('../data/standardized_log_exp_usage_matrix_by_v_gene_alpha_fmba.csv').drop(columns=['Unnamed: 0']).fillna(0)
    print(um.shape)
    process_all_files()
