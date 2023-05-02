import os
import time
from multiprocessing import Manager, Pool

import pandas as pd
from multipy.fwer import hochberg
from scipy.stats import fisher_exact

os.chdir('/home/ekvlasova/covid-t-cell-classifier')

hla_a = ['A*01', 'A*02', 'A*03', 'A*11', 'A*23', 'A*24', 'A*25', 'A*26', 'A*29', 'A*30', 'A*31',
         'A*32', 'A*33', 'A*66', 'A*68', 'A*69']

hla_b = ['B*07', 'B*08', 'B*13', 'B*14', 'B*15', 'B*18', 'B*27', 'B*35', 'B*37', 'B*38', 'B*39',
         'B*40', 'B*41', 'B*42', 'B*44', 'B*45', 'B*46', 'B*47', 'B*48', 'B*49', 'B*50',
         'B*51', 'B*52', 'B*53', 'B*54', 'B*55', 'B*56', 'B*57', 'B*58', 'B*73']
hla_a_to_results = Manager().dict()


def get_sum_clones_by_runs(runs):
    return run_to_number_of_clonotypes.loc[runs, :].number_of_clones.sum()


def get_sum_clone_occurences_with_hla_by_runs(runs, clone):
    return clonotype_matrix.loc[runs, :][clone].sum()


def process_one_clone_hla_a(params):
    clone, hla = params
    runs_with_hla = hla_to_runs_with_hla_a[hla]
    runs_without_hla = hla_to_runs_without_hla_a[hla]
    hla_clone = get_sum_clone_occurences_with_hla_by_runs(runs_with_hla, clone)
    hla_no_clone = get_sum_clones_by_runs(runs_with_hla) - hla_clone
    no_hla_clone = get_sum_clone_occurences_with_hla_by_runs(runs_without_hla, clone)
    no_hla_no_clone = get_sum_clones_by_runs(runs_without_hla) - no_hla_clone

    hla_a_to_results[clone] = fisher_exact([[hla_clone, hla_no_clone], [no_hla_clone, no_hla_no_clone]])[1]
    # if len(hla_a_to_results[hla]) % 1000 == 0:
    # print(f'dict size for {hla} is {len(hla_a_to_results)}')


def process_all_clones_hla(hla):
    start = time.time()
    print(f'Started with {hla}')
    with Pool(48) as p:
        p.map(process_one_clone_hla_a, [(x, hla) for x in clonotype_matrix.columns])
    print(f'Done with {hla} in {time.time() - start} seconds')


if __name__ == "__main__":
    run_to_number_of_clonotypes = pd.read_csv('data/run_to_number_of_clonotypes.csv').set_index('run')
    desc_fmba = pd.read_csv('/projects/fmba_covid/1_data_links/fmba_2021.txt', sep='\t')
    desc_fmba = desc_fmba[desc_fmba['file.name'].str.contains('TRB')].rename(columns={'file.name': 'run'})
    desc_hip = pd.read_csv('/projects/fmba_covid/1_data_links/hip_full_prep.txt', sep='\t')
    desc_hip = desc_hip.rename(columns={'file_name': 'run'})
    desc_adaptive = pd.read_csv('/projects/fmba_covid/1_data_links/adaptive_new.txt', sep='\t')
    desc_fmba = desc_fmba.dropna(subset=['HLA-A.1', 'HLA-A.2', 'HLA-B.1', 'HLA-B.2']).reset_index(drop=True)
    hla_list = []
    for hla_a1, hla_a2, hla_b1, hla_b2 in zip(desc_fmba['HLA-A.1'], desc_fmba['HLA-A.2'], desc_fmba['HLA-B.1'],
                                              desc_fmba['HLA-B.2']):
        ans = ''
        hla_a1 = hla_a1.split(':')[0]
        hla_a2 = hla_a2.split(':')[0]
        hla_b1 = hla_b1.split(':')[0]
        hla_b2 = hla_b2.split(':')[0]
        ans += hla_a1 + ','
        if hla_a1 != hla_a2:
            ans += hla_a2 + ','
        ans += hla_b1 + ','
        if hla_b1 != hla_b2:
            ans += hla_b2
        hla_list.append(ans)
    desc_fmba['hla'] = pd.Series(hla_list)
    um = pd.read_csv('data/standardized_log_exp_usage_matrix_joint_new.csv').drop(columns=['Unnamed: 0']).fillna(0)
    # um = um[(um.run.str.contains('HIP')) | (um.run.str.contains('Keck'))].reset_index(drop=True)
    all_hla_data = pd.concat([desc_hip[['run', 'hla']], desc_fmba[['run', 'hla']]]).dropna().merge(um)
    clonotype_matrix = pd.read_csv('data/covid_significant_clonotype_matrix_500k_top_0_mismatch.csv')
    clonotype_matrix = clonotype_matrix.drop(columns=['Unnamed: 0']).set_index('run')

    hla_to_runs_with_hla_a = {}
    hla_to_runs_without_hla_a = {}
    for hla in hla_b:
        hla_to_runs_with_hla_a[hla] = all_hla_data[all_hla_data.hla.str.contains(hla)].run
        hla_to_runs_without_hla_a[hla] = all_hla_data[~all_hla_data.hla.str.contains(hla)].run

    all_significant_clones = []
    for hla in hla_b:
        process_all_clones_hla(hla)

        clones = []
        pvals = []
        for clone, pval in hla_a_to_results.items():
            clones.append(clone)
            pvals.append(pval)
        significant_pvals = hochberg(pvals, alpha=0.05)
        significant_clones = []
        for pval, clone in zip(significant_pvals, clones):
            if pval:
                significant_clones.append(clone)
        all_significant_clones.append(
            pd.DataFrame(data={'clone': significant_clones, 'hla': [hla for _ in range(len(significant_clones))]}))

    pd.concat(all_significant_clones).to_csv('data/significant_hla_b_clones_500k_top_0mismatch.csv')