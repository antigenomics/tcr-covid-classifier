from multiprocessing import Manager, Pool

import numpy as np
import pandas as pd
from tqdm import tqdm
run_to_clone_number = Manager().dict()


def process_one_file(run):
    try:
        cur_data = pd.read_csv(f'{snakemake.input[1]}/{run}')
        # cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/hip_full/{run}', sep='\t')
        run_to_clone_number[run] = cur_data.shape[0]
        print(len(run_to_clone_number))
    except Exception as e:
        pass


def process_all_files():
    runs = list(um['run'])
    print(len(runs))
    with Pool(snakemake.threads) as p:
        p.map(process_one_file, runs)
    data_dict = {'run': [], 'number_of_clones': []}
    for x, y in run_to_clone_number.items():
        data_dict['run'].append(x)
        data_dict['number_of_clones'].append(y)
    pd.DataFrame(data=data_dict).to_csv(snakemake.output[0], index=False)


if __name__ == "__main__":
    um = pd.read_csv(snakemake.input[0]).drop(columns=['Unnamed: 0']).fillna(0)
    process_all_files()
