import pandas as pd

if __name__ == "__main__":
    desc = pd.read_csv('data/fmba_2021.txt', sep='\t')
    desc = desc[desc.COVID_status != 'unknown']
    desc.COVID_status = desc.COVID_status.apply(lambda x: x if x == 'COVID' else 'healthy')
    desc = desc[~(desc['HLA-A.1'].isna())]
    desc['run'] = desc['file.name']
    desc.to_csv('data/preprocessed_fmba_metadata.csv', index=False)
