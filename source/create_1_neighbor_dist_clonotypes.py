import os
from multiprocessing import Pool, Manager

import pandas as pd

os.chdir('/home/ekvlasova/covid-t-cell-classifier')

run_to_presence_of_clonotypes = Manager().dict()


def check_1_mismatch_neighboring_clones(clone1, clone2):
    if len(clone2) != len(clone1):
        return False
    n = len(clone1) // 2
    if n == 0:
        return True
    if clone1[:n] != clone2[:n] and clone1[n:] != clone2[n:]:
        return False
    return check_1_mismatch_neighboring_clones(clone1[:n], clone2[:n]) and check_1_mismatch_neighboring_clones(
        clone1[n:], clone2[n:])


def process_one_file(run):
    run, most_common_clonotypes, i = run
    try:
        cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/downsampled_new/{run}')
    except Exception:
        return
    cur_cdrs = set(cur_data['cdr3aa'])

    res = set()
    for clonotype in most_common_clonotypes['cdr3aa']:
        for cdr in cur_cdrs:
            if check_1_mismatch_neighboring_clones(clonotype, cdr):
                res.add(cdr)
    run_to_presence_of_clonotypes[run] = res

    if i % 300 == 0:
        print(f'processed {i} runs')


def process_all_files(save_path, most_common_clonotypes):
    desc = pd.read_csv('data/desc_fmba_not_nan_hla.csv')
    runs = [(x, most_common_clonotypes, i) for i, x in enumerate(desc['run'].tolist()) if 'TRB' in x]

    with Pool(80) as p:
        p.map(process_one_file, runs)

    all_clones = set()
    for clones in run_to_presence_of_clonotypes.values():
        all_clones.update(clones)

    pd.DataFrame.from_dict(data={'cdr3aa': pd.Series(list(all_clones))}).to_csv(save_path, index=False)


def find_neighbors_for_all_alleles():
    hla_keys = pd.read_csv('../data/hla_keys.csv')['0']
    print(hla_keys)
    for hla in list(hla_keys):
        print(f'Started processing {hla}')
        most_common_clonotypes = \
            pd.read_csv(f'../data/hla_top_clonotypes/most_used_500k_clonotypes_top_fmba_hla_{hla}.csv')[
                ['cdr3aa']].drop_duplicates().reset_index(drop=True)
        process_all_files(
            save_path=f'../data/hla_clonotype_matrix/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{hla}.csv',
            most_common_clonotypes=most_common_clonotypes)


if __name__ == "__main__":
    process_all_files(save_path=f'data/hla_covid_results/covid_hla_neighbors_500k_top_1_mismatch_hla_A*02.csv',
                      most_common_clonotypes=pd.read_csv(
                          f'data/hla_covid_results/covid_hla_clones_500k_top_1_mismatch_hla_A*02.csv').rename(
                          columns={'clone': 'cdr3aa'}))
