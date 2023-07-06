from multiprocessing import Manager, Pool

import numpy as np
import pandas as pd
from tqdm import tqdm
run_to_clone_number = Manager().dict()
run_to_read_number = Manager().dict()


def process_one_file(run):
    try:
        file = f'{snakemake.input[1]}/{run}'
        cur_data = pd.read_csv(file, sep=',' if file.split('.')[-1] == 'csv' else '\t')
        run_to_clone_number[run] = cur_data.shape[0]
        run_to_read_number[run] = cur_data['count'].sum()
    except Exception as e:
        pass


def process_all_files():
    if 'run' in um.columns:
        runs = list(um['run'])
    else:
        runs = list(um['file_name'])
    print(len(runs))
    with Pool(snakemake.threads) as p:
        p.map(process_one_file, runs)
    data_dict = {'run': [], 'number_of_clones': []}
    if snakemake.params.make_read_count_col:
        data_dict['number_of_reads'] =  []
    for x, y in run_to_clone_number.items():
        data_dict['run'].append(x)
        data_dict['number_of_clones'].append(y)
        if snakemake.params.make_read_count_col:
            data_dict['number_of_reads'].append(run_to_read_number[x])
    pd.DataFrame(data=data_dict).to_csv(snakemake.output[0], index=False)


if __name__ == "__main__":
    file = snakemake.input[0]
    um = pd.read_csv(file, sep=',' if file.split('.')[-1] == 'csv' else '\t')
    print(um)
    process_all_files()
