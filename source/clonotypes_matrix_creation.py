import os
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
            occurences.add(clone_to_search)
    return occurences


def check(clone1, clone2, mismatch_max=0):
    ans = 0
    for c1, c2 in zip(clone1, clone2):
        if c1 != c2:
            ans += 1
    return ans <= mismatch_max


def process_one_file(run):
    run, mismatch_max, most_common_clonotypes, raw_data_folder, index = run
    try:
        cur_data = pd.read_csv(f'{raw_data_folder}/{run}')
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

    except Exception as e:
        print(f'error processing {raw_data_folder}/{run}', e)
        return

    if index % 200 == 0:
        print(f'processed {index} runs')


def process_all_files(save_path, most_common_clonotypes, um, mismatch_max=0, raw_data_folder='downsampled_new'):
    run_to_presence_of_clonotypes['cdr3aa'] = most_common_clonotypes['cdr3aa']
    runs = [(x, mismatch_max, most_common_clonotypes, raw_data_folder, i) for i, x in enumerate(um['run'].tolist())]
    print('Start processing')
    with Pool(snakemake.threads) as p:
        p.map(process_one_file, runs)#, chunksize=10)
        # p.close()
        # p.terminate()
        # p.join()
    print(run_to_presence_of_clonotypes)
    data = {x: y for x, y in run_to_presence_of_clonotypes.items()}
    pd.DataFrame.from_dict(data=data).to_csv(save_path, index=False)


def create_matrix_for_allele_data(um_path, top_clonotypes_path, clone_matrices_path, raw_data_folder, hla_keys=None):
    um = pd.read_csv(um_path).drop(columns=['Unnamed: 0']).fillna(0)
    os.mkdir(clone_matrices_path)
    if hla_keys is None:
        hla_keys = pd.read_csv('data/hla_keys.csv')['0']
    for hla in list(hla_keys):
        print(f'Started processing {hla}')
        most_common_clonotypes = \
            pd.read_csv(f'{top_clonotypes_path}/most_used_clonotypes_fmba_hla_{hla}.csv')[
                ['cdr3aa']].drop_duplicates().reset_index(drop=True)
        process_all_files(
            save_path=f'{clone_matrices_path}/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{hla}.csv',
            most_common_clonotypes=most_common_clonotypes,
            um=um,
            mismatch_max=1,
            raw_data_folder=raw_data_folder)


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
        if snakemake.params.platform == 'adaptive':
            um = pd.read_csv(snakemake.input[0]).drop(columns=['Unnamed: 0']).fillna(0)
            process_all_files(save_path=snakemake.output[0],
                              most_common_clonotypes=pd.read_csv(snakemake.input[2]),
                              um=um,
                              mismatch_max=1,
                              raw_data_folder=snakemake.input[1],
                              )
        if snakemake.params.platform == 'fmba-allele':
            create_matrix_for_allele_data(um_path=snakemake.input[0],
                                          top_clonotypes_path=snakemake.input[1],
                                          clone_matrices_path=snakemake.output[0],
                                          raw_data_folder=snakemake.input[2],
                                          hla_keys=snakemake.params.hla_to_consider)

