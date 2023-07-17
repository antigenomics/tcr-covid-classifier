from collections import defaultdict
import pandas as pd

desc_fmba = pd.read_csv(snakemake.input[0])
desc_fmba = desc_fmba[desc_fmba['run'].str.contains('TRB')].rename(columns={'file.name': 'run'})

hla_to_patients = defaultdict(set)
for hla_variant in ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPB1', 'HLA-DQB1', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB3', 'HLA-DRB5']:
    if hla_variant + '.1' in desc_fmba.columns:
        cur_desc = desc_fmba.dropna(subset=[hla_variant + '.1'])
        for run, hla_res in zip(cur_desc['run'], cur_desc[hla_variant + '.1']):
            hla_to_patients[hla_res.split(':')[0]].add(run)
    if hla_variant + '.2' in desc_fmba.columns:
        cur_desc = desc_fmba.dropna(subset=[hla_variant + '.2'])
        for run, hla_res in zip(cur_desc['run'], cur_desc[hla_variant + '.2']):
            hla_to_patients[hla_res.split(':')[0]].add(run)

hla_keys = list(hla_to_patients.keys())

pd.Series(hla_to_patients.keys()).to_csv(snakemake.output[0])

import os

os.mkdir(snakemake.output[1])

for hla in hla_to_patients:
    current_subdata = desc_fmba[desc_fmba.run.isin(hla_to_patients[hla])].reset_index(drop=True)
    current_subdata['covid'] = current_subdata.COVID_status.apply(lambda x: 'covid' if x == 'COVID' else 'healthy')
    current_subdata['is_test_run'] = current_subdata.folder.str.lower().str.contains('novaseq6')
    current_subdata.to_csv(f'{snakemake.output[1]}/fmba_desc_hla_{hla}.csv')
