import numpy as np
import pandas as pd
from tqdm import tqdm

hla_a = ['A*01',
         'A*02',
         'A*03',
         'A*11',
         'A*23',
         'A*24',
         'A*25',
         'A*26',
         'A*29',
         'A*30',
         'A*31',
         'A*32',
         'A*33',
         'A*66',
         'A*68',
         'A*69']

hla_b = ['B*07',
         'B*08',
         'B*13',
         'B*14',
         'B*15',
         'B*18',
         'B*27',
         'B*35',
         'B*37',
         'B*38',
         'B*39',
         'B*40',
         'B*41',
         'B*42',
         'B*44',
         'B*45',
         'B*46',
         'B*47',
         'B*48',
         'B*49',
         'B*50',
         'B*51',
         'B*52',
         'B*53',
         'B*54',
         'B*55',
         'B*56',
         'B*57',
         'B*58',
         'B*73']


def process_all_files(usage_matrix_path):
    um = pd.read_csv(usage_matrix_path).drop(columns=['Unnamed: 0']).fillna(0).merge(
        pd.read_csv('/projects/fmba_covid/1_data_links/fmba_2021.txt', sep='\t').rename(
            columns={'file.name': 'run'})[['run', 'HLA-A.1', 'HLA-A.2', 'HLA-B.1', 'HLA-B.2']]
    )
    datasets_to_concat = []
    for run, covid, hla_a1, hla_a2, hla_b1, hla_b2 in tqdm(
            zip(um['run'], um['covid'], um['HLA-A.1'], um['HLA-A.2'], um['HLA-B.1'], um['HLA-B.2']), total=len(um),
            desc='Extracting top clonotypes'):
        try:
            cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/downsampled_new/{run}', nrows=15000)
            cur_data['covid'] = 1 if covid == 'covid' else None
            cur_data['healthy'] = 1 if covid == 'healthy' else None
            for hla in hla_a:
                cur_data[hla] = 1 if hla in hla_a1 or hla in hla_a2 else None
            for hla in hla_b:
                cur_data[hla] = 1 if hla in hla_b1 or hla in hla_b2 else None
            datasets_to_concat.append(cur_data)
        except Exception as e:
            pass
    full_data = pd.concat(datasets_to_concat)
    top = full_data.groupby(['cdr3aa', 'v', 'j'], as_index=False).count()
    top = top.sort_values(by=['count'], ascending=False).head(10000)
    top.to_csv(f'../data/most_used_clonotypes_covid_fmba_{10000}.csv', index=False)


if __name__ == "__main__":
    process_all_files('../data/standardized_log_exp_usage_matrix_functional.csv')
