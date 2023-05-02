import os

os.chdir('/home/ekvlasova/covid-t-cell-classifier')

import pandas as pd


def create_allele_column(desc, allele):
    runs_with_allele = pd.read_csv(f'data/hla_desc/fmba_desc_hla_{allele}.csv').run
    runs_with_allele = set(runs_with_allele.apply(lambda x: x.split('.')[0]))
    desc[allele] = desc.run.apply(lambda x: 'present' if x in runs_with_allele else 'absent')
    return desc


def evaluate_desc(desc_path, hla_keys):
    desc = pd.read_csv(desc_path).drop(columns=['Unnamed: 0'])
    desc['run'] = desc['run'].apply(lambda x: x.split('.')[0])
    for hla in hla_keys:
        desc = create_allele_column(desc, hla)
    return desc


if __name__ == "__main__":
    hla_keys = pd.read_csv('data/hla_keys.csv')['0']
    desc = evaluate_desc(desc_path=f'data/desc_fmba_not_nan_hla.csv',
                         hla_keys=hla_keys)
    desc.to_csv('data/desc_fmba_hla_bool.csv')
