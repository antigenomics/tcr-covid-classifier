import os

os.chdir('/home/ekvlasova/covid-t-cell-classifier')

import pandas as pd
from tqdm import tqdm
from multiprocessing import Manager, Pool
from utils.data_utils import prepare_clonotype_matrix

clone_to_df = Manager().dict()


def save_current_data(save_path):
    datasets_to_concat = []
    for clone, found_data in clone_to_df.items():
        datasets_to_concat.append(found_data)

    pd.concat(datasets_to_concat).to_csv(save_path)


def find_associated_alpha_clones_for_beta_clone(beta_clone, beta_cm, alpha_cm, save_path, alpha=0.95):
    runs_with_beta_clone = beta_cm[beta_cm[beta_clone] != 0].run
    runs_without_beta_clone = beta_cm[beta_cm[beta_clone] == 0].run
    alpha_cm_for_beta_clone = alpha_cm[alpha_cm.run.isin(runs_with_beta_clone)].drop(columns=['run'])
    usages = alpha_cm_for_beta_clone.sum(axis=0) / len(runs_with_beta_clone)
    usages = pd.DataFrame(usages, columns=['usage']).reset_index().rename(columns={'index': 'alpha_clone'})
    usages = usages[usages.usage > alpha]
    usages['beta_clone'] = beta_clone
    clone_to_df[beta_clone] = usages
    clone_to_df_size = len(clone_to_df)
    print(f'Processed {beta_clone}, iteration {clone_to_df_size}, found {len(usages)} associated alphas')
    if clone_to_df_size % 1000 == 0:
        save_current_data(save_path)
    return usages


def mapping_function(args):
    beta_clone, beta_cm, alpha_cm, save_path, alpha = args
    return find_associated_alpha_clones_for_beta_clone(beta_clone, beta_cm, alpha_cm, alpha)


def create_pairs_data(alpha_clone_matrix_path, beta_clone_matrix_path, save_path, alpha=0.95):
    datasets_to_concat = []
    alpha_cm = prepare_clonotype_matrix(alpha_clone_matrix_path, make_bool_features=True)
    beta_cm = prepare_clonotype_matrix(beta_clone_matrix_path, make_bool_features=True)
    print(f'Started searching for associations')
    arguments = []
    for i, clone in tqdm(enumerate(beta_cm.columns)):
        if clone == 'run':
            continue
        arguments.append([clone, beta_cm, alpha_cm, save_path, alpha])

    with Pool(40) as p:
        p.map(mapping_function, arguments)

    for clone, found_data in clone_to_df.items():
        datasets_to_concat.append(found_data)

    pd.concat(datasets_to_concat).to_csv(save_path)


if __name__ == "__main__":
    create_pairs_data(alpha_clone_matrix_path='data/alpha/clonotype_matrix_fmba_alpha_top_500k_0_mismatch.csv',
                      beta_clone_matrix_path='data/hla_sign_clone_matrix/covid_clonotype_matrix_500k_top_1_mismatch_hla_all_fisher.csv',
                      save_path='data/alpha_beta_pairs/paired_clones_alpha_mis_0_beta_mis_1_alpha_095_covid.csv',
                      alpha=0.95)
