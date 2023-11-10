import pandas as pd
from tqdm import tqdm
import os


def prepare_run_column(df):
    df['run'] = df['run'].apply(lambda x: x.split('.')[0])
    return df


def process_all_files(usage_matrix_path, save_path, method='top', count_of_clones=500000,
                      raw_data_folder='downsampled_new', runs_to_process=None):
    um = pd.read_csv(usage_matrix_path).drop(columns=['Unnamed: 0', 'covid']).fillna(0)
    datasets_to_concat = []

    for run in tqdm(um['run'], total=len(um), desc='Extracting top clonotypes'):
        if runs_to_process is not None and run.split('.')[0] not in runs_to_process:
            continue
        try:
            cur_data = pd.read_csv(f'{raw_data_folder}/{run}')
            cur_data = cur_data[['cdr3aa']].drop_duplicates()
            cur_data['count'] = 1
            datasets_to_concat.append(cur_data)
        except Exception as e:
            continue
    full_data = pd.concat(datasets_to_concat)
    full_data = full_data[full_data.cdr3aa.str.isalpha()]
    print(len(full_data))
    top = full_data.groupby(['cdr3aa'], as_index=False).count()

    if method == 'top':
        top = top[top['count'] < 0.4 * len(um)].sort_values(by=['count'], ascending=False).head(count_of_clones)

    elif method == 'random-roulette':
        top = top.sample(n=count_of_clones, random_state=42, weights='count')

    elif method == 'random-uniform':
        top = top.sample(n=count_of_clones, random_state=42)

    elif method == 'unique-occurence':
        top = top[top['count'] >= count_of_clones]
    top.to_csv(save_path, index=False)


def clonotype_extraction_for_allele(path_to_hla_desc, um_path, top_clonotypes_path, raw_data_folder, hla_keys=None):
    if hla_keys is None:
        hla_keys = pd.read_csv('data/hla_keys.csv')['0']
    print(hla_keys)
    os.mkdir(top_clonotypes_path)
    for hla in hla_keys:
        print(f'Started processing {hla}')
        hla_desc = pd.read_csv(f'{path_to_hla_desc}/fmba_desc_hla_{hla}.csv')
        runs_to_process = list(prepare_run_column(hla_desc).run)
        process_all_files(usage_matrix_path=um_path,
                          save_path=f'{top_clonotypes_path}/most_used_clonotypes_fmba_hla_{hla}.csv',
                          method='top',
                          count_of_clones=500000,
                          raw_data_folder=raw_data_folder,
                          runs_to_process=runs_to_process)


def clonotype_extraction_for_covid_healthy():
    fmba_desc = prepare_run_column(pd.read_csv('../data/desc_fmba_not_nan_hla.csv').drop(columns=['Unnamed: 0']))
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=f'../data/covid_top_clonotypes/most_used_500k_clonotypes_top_fmba_covid.csv',
                      method='top',
                      count_of_clones=500000,
                      runs_to_process=list(fmba_desc[fmba_desc.COVID_status == 'COVID'].run))
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=f'../data/covid_top_clonotypes/most_used_500k_clonotypes_top_fmba_healthy.csv',
                      method='top',
                      count_of_clones=500000,
                      runs_to_process=list(fmba_desc[fmba_desc.COVID_status != 'COVID'].run))


def clonotypes_extraction_procedure_for_adaptive():
    desc = prepare_run_column(pd.read_csv('../data/standardized_log_exp_usage_matrix_joint_new.csv'))
    desc = desc[desc.platform == 'adaptive']
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=f'../data/adaptive_clone_results/most_used_500k_clonotypes_top_adaptive.csv',
                      method='top',
                      count_of_clones=500000,
                      runs_to_process=list(desc.run))


def clonotypes_extraction_procedure_for_adaptive_batch():
    desc = pd.read_csv('../data/desc_hip_bool.csv')
    desc = desc[desc.platform == 'adaptive']
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=f'../data/adaptive_clone_results/most_used_100k_clonotypes_top_healthy.csv',
                      method='top',
                      count_of_clones=100000,
                      runs_to_process=list(desc[desc.is_keck == 'yes'].run))
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=f'../data/adaptive_clone_results/most_used_100k_clonotypes_top_ill.csv',
                      method='top',
                      count_of_clones=100000,
                      runs_to_process=list(desc[desc.is_keck != 'yes'].run))


def clonotypes_extraction_procedure_for_fmba(um_path, save_path, n_clones, resampled_samples_path, method):
    desc = prepare_run_column(pd.read_csv(um_path))
    process_all_files(usage_matrix_path=um_path,
                      save_path=save_path,
                      method=method,
                      count_of_clones=n_clones,
                      raw_data_folder=resampled_samples_path,
                      runs_to_process=list(desc.run))


def clonotype_extraction_procedure_by_batch_list(batch_list, save_path):
    desc = prepare_run_column(pd.read_csv('../data/standardized_log_exp_usage_matrix_joint_new.csv'))
    desc = desc[desc.project.isin(batch_list)]
    process_all_files(usage_matrix_path='../data/standardized_log_exp_usage_matrix_joint_new.csv',
                      save_path=save_path,
                      method='top',
                      count_of_clones=500000,
                      runs_to_process=list(desc.run))


if __name__ == "__main__":
    if 'snakemake' in globals():
        if snakemake.params.platform == 'fmba':
            clonotypes_extraction_procedure_for_fmba(um_path=snakemake.input[0],
                                                     save_path=snakemake.output[0],
                                                     method=snakemake.params.sampling_method,
                                                     n_clones=snakemake.params.n_clones,
                                                     resampled_samples_path=snakemake.input[1])
        if snakemake.params.platform == 'adaptive':
            clonotypes_extraction_procedure_for_fmba(um_path=snakemake.input[0],
                                                     save_path=snakemake.output[0],
                                                     method=snakemake.params.sampling_method,
                                                     n_clones=snakemake.params.n_clones,
                                                     resampled_samples_path=snakemake.input[1])
        if snakemake.params.platform == 'fmba-allele':
            clonotype_extraction_for_allele(path_to_hla_desc=snakemake.input[0],
                                            um_path=snakemake.input[1],
                                            raw_data_folder=snakemake.input[2],
                                            top_clonotypes_path=snakemake.output[0],
                                            hla_keys=snakemake.params.hla_to_consider)
