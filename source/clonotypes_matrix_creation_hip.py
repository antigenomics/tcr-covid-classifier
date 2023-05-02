import time
from collections import defaultdict, Counter
from multiprocessing import Pool, Manager

import numpy as np
import pandas as pd
from tqdm import tqdm

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


def get_count_occurences(data, clone):
    def check_identity(clone1, clone2):
        if len(clone1) != len(clone2):
            return False
        has_mismatches = False
        for c1, c2 in zip(clone1, clone2):
            if c1 != c2 and has_mismatches:
                return False
            elif c1 != c2:
                has_mismatches = True
        return True
    return len(data[data.cdr3aa.apply(lambda x: check_identity(x, clone))])


def process_one_file(run):
    start = time.time()
    try:
        # cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/downsampled_new/{run}')
        cur_data = pd.read_csv(f'/projects/fmba_covid/1_data_links/hip_full/{run}', sep='\t')
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
                if check(clone, cdr3_mismatch_clone, mismatch_max=0):
                    sum_occurences += cur_cdrs[cdr3_mismatch_clone]

            res.append(sum_occurences)
        else:
            res.append(0)
    run_to_presence_of_clonotypes[run] = pd.Series(res)
    print(f'finished {run} in {time.time() - start} seconds')


def process_all_files():
    run_to_presence_of_clonotypes['cdr3aa'] = most_common_clonotypes['cdr3aa']
    runs = um['run'].tolist()
    with Pool(48) as p:
        p.map(process_one_file, runs)
    data = {x: y for x, y in run_to_presence_of_clonotypes.items()}
    pd.DataFrame.from_dict(data=data).to_csv(f'../data/clonotype_matrix_hip_10k_0_mismatch_roulette.csv', index=False)


if __name__ == "__main__":
    um = pd.read_csv('../data/standardized_log_exp_usage_matrix_joint_new.csv').drop(columns=['Unnamed: 0']).fillna(0)
    um = um[(um.run.str.contains('HIP')) | (um.run.str.contains('Keck'))]
    most_common_clonotypes = pd.read_csv(f'../data/most_used_clonotypes_hip_10000_wo_nt_roulette.csv')[['cdr3aa']]
    process_all_files()
