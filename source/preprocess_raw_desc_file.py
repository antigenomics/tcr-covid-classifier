import pandas as pd

if __name__ == "__main__":
    desc = pd.read_csv(snakemake.input[0], sep='\t')
    desc = desc[desc.COVID_status != 'unknown']
    desc.COVID_status = desc.COVID_status.apply(lambda x: x if x == 'COVID' else 'healthy')
    desc = desc[~(desc['HLA-A.1'].isna())]
    desc['run'] = desc['file.name']
    # desc = desc[(desc['HLA-A.1'].str.contains(r'A\*02:*')) | (desc['HLA-A.2'].str.contains(r'A\*02:*'))]
    print(desc.shape)
    desc.to_csv(snakemake.output[0], index=False)
