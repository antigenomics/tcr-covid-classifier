import pandas as pd


def prepare_run_column(df):
    df['run'] = df['run'].apply(lambda x: x.split('.')[0])
    return df


def prepare_clonotype_matrix(clonotype_matrix_path, make_bool_features=False):
    print(f'Preparing {clonotype_matrix_path}')
    df = pd.read_csv(clonotype_matrix_path)
    if 'Unnamed: 0' in df.columns:
        df = df.drop(columns=['Unnamed: 0'])
    if 'cdr3aa' in df.columns:
        df = df.set_index('cdr3aa').T.reset_index().rename(columns={'index': 'run'})
    df = prepare_run_column(df)
    run_col = df.run
    if make_bool_features:
        df = df.drop(columns=['run'])
        df = df.astype(bool)
        df['run'] = run_col
    return df


def get_anomaly_clones(clonotype_matrix_path, desc_path=None, critical_percent=0.4):
    cm = prepare_clonotype_matrix(clonotype_matrix_path, make_bool_features=True)
    if desc_path is not None:
        cm = cm.merge(prepare_clonotype_matrix(desc_path)[['run']])
    results = pd.DataFrame(cm.drop(columns=['run']).sum(axis=0)).reset_index().rename(columns={'index': 'clone', 0: 'count'})
    return results[results['count'] > cm.shape[0] * critical_percent]
