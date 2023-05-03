import time
from collections import defaultdict, Counter
from multiprocessing import Pool, Manager

import pandas as pd

run_to_presence_of_clonotypes = Manager().dict()

def check_mismatch_clone(cur_clone, mismatch_clones):
    occurences = set()
    for i in range(len(cur_clone)):
        clone_to_search = cur_clone[: i] + 'X' + cur_clone[i + 1:]
        cur_clones_set = mismatch_clones[i]
        if clone_to_search in cur_clones_set:
            # print(cur_clone, clone_to_search, mismatch_clones[i])
            occurences.add(clone_to_search)
    return occurences


def check(clone1, clone2, mismatch_max=0):
    ans = 0
    for c1, c2 in zip(clone1, clone2):
        if c1 != c2:
            ans += 1
    return ans <= mismatch_max


def process_one_file(run):
    run, mismatch_max, most_common_clonotypes, raw_data_folder, i = run
    start = time.time()
    try:
        cur_data = pd.read_csv(f'{raw_data_folder}/{run}')
    except Exception as e:
        return
    res = []
    cur_cdrs = Counter(cur_data['cdr3aa'])
    length_to_clones = defaultdict(set)
    for cdr in cur_cdrs:
        length_to_clones[len(cdr)].add(cdr)

    length_to_mismatch_clones = {}
    mismatch_clone_to_cdr3aa = defaultdict(set)
    for length, cdr_set in length_to_clones.items():
        length_to_mismatch_clones[length] = defaultdict(set)
        for clone in cdr_set:
            if not clone.isalpha():
                continue
            for i in range(len(clone)):
                mismatch_clone = clone[:i] + 'X' + clone[i + 1:]
                length_to_mismatch_clones[length][i].add(mismatch_clone)
                mismatch_clone_to_cdr3aa[mismatch_clone].add(clone)

    for clone in most_common_clonotypes['cdr3aa']:
        if len(clone) in length_to_mismatch_clones:
            mismatch_clones = check_mismatch_clone(clone, length_to_mismatch_clones[len(clone)])
            cdr3aa_found_clones = set()
            for mismatch_clone in mismatch_clones:
                cdr3aa_found_clones.update(mismatch_clone_to_cdr3aa[mismatch_clone])
            sum_occurences = 0
            for cdr3_mismatch_clone in cdr3aa_found_clones:
                if check(clone, cdr3_mismatch_clone, mismatch_max=mismatch_max):
                    sum_occurences += cur_cdrs[cdr3_mismatch_clone]

            res.append(sum_occurences)
        else:
            res.append(0)
    run_to_presence_of_clonotypes[run] = pd.Series(res)

    # print(f'finished {run} in {time.time() - start} seconds')
    if i % 500 == 0:
        print(f'processed {i} runs')


def process_all_files(save_path, most_common_clonotypes, um, mismatch_max=0, raw_data_folder='downsampled_new'):
    run_to_presence_of_clonotypes['cdr3aa'] = most_common_clonotypes['cdr3aa']
    runs = [(x, mismatch_max, most_common_clonotypes, raw_data_folder, i) for i, x in enumerate(um['run'].tolist())]
    # runs = [(x, mismatch_max) for x in um['run'].tolist() if 'Keck' in x or 'HIP' in x]
    with Pool(snakemake.threads) as p:
        p.map(process_one_file, runs)
    data = {x: y for x, y in run_to_presence_of_clonotypes.items()}
    # pd.DataFrame.from_dict(data=data).to_csv(f'../data/clonotype_matrix_50k_{mismatch_max}_mismatch_top.csv', index=False)
    pd.DataFrame.from_dict(data=data).to_csv(save_path, index=False)


def create_matrix_for_allele_data():
    um = pd.read_csv('../data/standardized_log_exp_usage_matrix_joint_new.csv').drop(columns=['Unnamed: 0']).fillna(0)
    hla_keys = pd.read_csv('../data/hla_keys.csv')['0']
    for hla in list(hla_keys):
        print(f'Started processing {hla}')
        most_common_clonotypes = \
            pd.read_csv(f'../data/hla_top_clonotypes/most_used_500k_clonotypes_top_fmba_hla_{hla}.csv')[
                ['cdr3aa']].drop_duplicates().reset_index(drop=True)
        process_all_files(
            save_path=f'../data/hla_clonotype_matrix/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{hla}.csv',
            most_common_clonotypes=most_common_clonotypes,
            um=um,
            mismatch_max=1)


def clonotype_matrix_for_beta_chain():
    um = pd.read_csv('../data/standardized_log_exp_usage_matrix_joint_new.csv').drop(columns=['Unnamed: 0']).fillna(0)
    um = um[um.platform == 'adaptive']
    process_all_files(save_path='../data/clonotype_matrix_adaptive_top_500k_0_mismatch.csv',
                      most_common_clonotypes=pd.read_csv(f'../data/adaptive_clone_results/most_used_500k_clonotypes_top_adaptive.csv'),
                      um=um,
                      mismatch_max=0)


def clonotype_matrix_for_alpha_chain():
    um = pd.read_csv('../data/standardized_log_exp_usage_matrix_by_v_gene_alpha_fmba.csv').drop(
        columns=['Unnamed: 0']).fillna(0)
    process_all_files(save_path='../data/alpha/clonotype_matrix_fmba_alpha_top_500k_0_mismatch.csv',
                      most_common_clonotypes=pd.read_csv(
                          f'../data/alpha/most_used_500k_clonotypes_top.csv'),
                      raw_data_folder='downsampled_alpha',
                      um=um,
                      mismatch_max=0)


def clonotype_matrix_for_projects(projects_list, most_common_clones_path, save_path, mismatch=0):
    um = pd.read_csv('../data/standardized_log_exp_usage_matrix_joint_new.csv').drop(columns=['Unnamed: 0']).fillna(0)
    um = um[um.project.isin(projects_list)]
    process_all_files(save_path=save_path,
                          most_common_clonotypes=pd.read_csv(most_common_clones_path),
                          um=um,
                          mismatch_max=mismatch)


if __name__ == "__main__":
    if 'snakemake' in globals():
        if snakemake.params.platform == 'fmba':
            um = pd.read_csv(snakemake.input[0]).drop(columns=['Unnamed: 0']).fillna(0)
            process_all_files(save_path=snakemake.output[0],
                              most_common_clonotypes=pd.read_csv(snakemake.input[2]),
                              um=um,
                              mismatch_max=1,
                              raw_data_folder=snakemake.input[1],
                              )