import pandas as pd


if __name__ == "__main__":
    df = pd.read_csv('/projects/fmba_covid/1_data_links/hip_full.txt', sep='\t')
    df['dataset'] = df['sample_id'].apply(lambda x: 'HIP' if x.startswith('HIP') else 'KECK')
    df['file_name'] = df['file_name'].apply(lambda x: x.split('/')[-1])
    df['covid'] = 'healthy'
    df.to_csv('/projects/fmba_covid/1_data_links/hip_full_prep.txt', sep='\t', index=False)
