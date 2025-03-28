import pandas as pd
from tqdm import tqdm


def generate_v_gene_usage_vector(sample, functional_seqs_type='all', gene_type='v'):
    sample[f'{gene_type}_gene'] = sample[gene_type].apply(lambda x: x.split(',')[0])
    working_data = sample[[f'{gene_type}_gene', 'cdr3nt', 'cdr3aa']]
    if functional_seqs_type == 'functional':
        working_data = working_data[working_data['cdr3aa'].str.isalpha()]
    elif functional_seqs_type == 'nonfunctional':
        working_data = working_data[~working_data['cdr3aa'].str.isalpha()]
    v_gene_usage_df = working_data.sort_values(by=f'{gene_type}_gene').groupby(f'{gene_type}_gene',
                                                                               as_index=False).count()
    gene_to_usage_dict = {}
    for k, v in zip(v_gene_usage_df[f'{gene_type}_gene'], v_gene_usage_df['cdr3nt']):
        if k.startswith(starting_pattern):
            gene_to_usage_dict[k] = v
    return gene_to_usage_dict


def process_all_files(all_runs_path, names_mapping, raw_data_path, dataset='fmba', get_extra_info=False,
                      chain_to_read='TRB', functional_seqs_type='all', gene_type='v'):
    all_v_genes = set()
    runs = pd.read_csv(all_runs_path)
    sep = '\t'
    if dataset == 'fmba':
        sep = ','
        runs = runs[(runs[names_mapping['file_name']].str.contains(chain_to_read))]
    elif dataset == 'adaptive_new':
        runs[names_mapping['file_name']] = runs[names_mapping['file_name']]
    elif dataset == 'hip_full':
        runs[names_mapping['file_name']] = runs[names_mapping['file_name']].apply(lambda x: x.split('/')[-1])
    repertoires_list = []
    i = 0
    for run, project in tqdm(zip(runs[names_mapping['file_name']], runs[names_mapping['dataset']]), total=len(runs),
                             desc='Extracting v gene usage vectors'):
        data = pd.read_csv(f'{raw_data_path}/{dataset}/{run}', sep=sep)
        run_repertoire = generate_v_gene_usage_vector(data, functional_seqs_type=functional_seqs_type,
                                                      gene_type=gene_type)
        all_v_genes.update(list(run_repertoire.keys()))
        repertoires_list.append(run_repertoire)

    full_usage_dict = {}
    for v_gene in tqdm(all_v_genes, desc='Making up usage dataframe'):
        v_gene_vector = []
        for repertoire_dict in repertoires_list:
            if v_gene in repertoire_dict:
                v_gene_vector.append(repertoire_dict[v_gene])
            else:
                v_gene_vector.append(0)
        full_usage_dict[v_gene] = v_gene_vector

    full_usage_dict['run'] = runs[names_mapping['file_name']]
    full_usage_dict['project'] = runs[names_mapping['dataset']]
    runs[names_mapping['covid']] = runs[names_mapping['covid']].fillna('healthy')
    full_usage_dict['covid'] = runs[names_mapping['covid']].apply(
        lambda x: 'covid' if 'covid' in x.lower() else 'healthy')
    if get_extra_info:
        full_usage_dict['hla'] = runs['HLA-A.2']
        full_usage_dict['number_of_clonotypes'] = runs['clonotypes']

    usage_df = pd.DataFrame.from_dict(full_usage_dict)
    print(usage_df.columns)
    return usage_df
    # usage_df.to_csv(output_file)


def run_joint_matrix_creation():
    fmba = process_all_files(
        all_runs_path='/projects/fmba_covid/1_data_links/fmba_2021.txt',
        names_mapping={'file_name': 'file.name',
                       'covid': 'COVID_status',
                       'dataset': 'folder'},
        dataset='fmba_2021',
        # get_extra_info=True,
        get_extra_info=False,
        chain_to_read='TRB'
    )
    # ).to_csv('../data/usage_matrix.csv')
    adaptive = process_all_files(
        all_runs_path='/projects/fmba_covid/1_data_links/adaptive_new.txt',
        names_mapping={'file_name': 'sample_name',
                       'covid': 'Virus Diseases',
                       'dataset': 'Dataset'},
        dataset='adaptive_new',
        chain_to_read='TRB'
    )
    hip = process_all_files(
        all_runs_path='/projects/fmba_covid/1_data_links/hip_full_prep.txt',
        names_mapping={'file_name': 'file_name',
                       'covid': 'covid',
                       'dataset': 'dataset'},
        dataset='hip_full',
        chain_to_read='TRB'
    )
    pd.concat([fmba, adaptive, hip]).to_csv('../data/usage_matrix_joint_new.csv')


def adaptive_hip_matrix_creation(hip_desc_path, adaptive_desc_path, raw_path, save_path, gene_type):
    adaptive = process_all_files(
        all_runs_path=adaptive_desc_path,
        names_mapping={'file_name': 'file_name',
                       'covid': 'covid',
                       'dataset': 'Dataset'},
        dataset='adaptive_new',
        chain_to_read='TRB',
        raw_data_path=raw_path,
        gene_type=gene_type
    )
    hip = process_all_files(
        all_runs_path=hip_desc_path,
        names_mapping={'file_name': 'file_name',
                       'covid': 'covid',
                       'dataset': 'dataset'},
        dataset='hip_full',
        chain_to_read='TRB',
        raw_data_path=raw_path,
        gene_type=gene_type
    )
    bad_coverage_files = pd.read_csv('data/bad_quality_adaptive_runs.csv').run
    final_data = pd.concat([adaptive, hip])
    final_data[~final_data.run.isin(bad_coverage_files)].to_csv(save_path)


def run_alpha_chain_fmba_matrix_creation(desc_path, output_path, raw_data_path, gene_type):
    process_all_files(
        all_runs_path=desc_path,
        names_mapping={'file_name': 'run',
                       'covid': 'COVID_status',
                       'dataset': 'folder'},
        raw_data_path=raw_data_path,
        dataset='fmba',
        get_extra_info=False,
        chain_to_read='TRA',
        gene_type=gene_type
    ).to_csv(output_path)


def run_beta_chain_fmba_matrix_creation(desc_path, output_path, raw_data_path, gene_type):
    process_all_files(
        all_runs_path=desc_path,
        names_mapping={'file_name': 'run',
                       'covid': 'COVID_status',
                       'dataset': 'folder'},
        raw_data_path=raw_data_path,
        dataset='fmba',
        get_extra_info=False,
        chain_to_read='TRB',
        gene_type=gene_type
    ).to_csv(output_path)


def run_beta_chain_adaptive_matrix_creation():
    process_all_files(
        all_runs_path='/projects/fmba_covid/1_data_links/adaptive_new.txt',
        names_mapping={'file_name': 'sample_name',
                       'covid': 'Virus Diseases',
                       'dataset': 'Dataset'},
        dataset='adaptive_new',
        chain_to_read='TRB'
    ).to_csv('../data/usage_matrix_adaptive.csv')


def run_beta_chain_hip_matrix_creation():
    process_all_files(
        all_runs_path='/projects/fmba_covid/1_data_links/hip_full_prep.txt',
        names_mapping={'file_name': 'file_name',
                       'covid': 'covid',
                       'dataset': 'dataset'},
        dataset='hip_full',
        chain_to_read='TRB'
    ).to_csv('../data/usage_matrix_hip.csv')


def run_beta_chain_fmba_matrix_creation_functional_nonfunctional():
    process_all_files(
        all_runs_path='/projects/fmba_covid/1_data_links/fmba_2021.txt',
        names_mapping={'file_name': 'file.name',
                       'covid': 'COVID_status',
                       'dataset': 'folder'},
        dataset='fmba_2021',
        get_extra_info=False,
        chain_to_read='TRB',
        functional_seqs_type='functional'
    ).to_csv('../data/usage_matrix_functional.csv')
    process_all_files(
        all_runs_path='/projects/fmba_covid/1_data_links/fmba_2021.txt',
        names_mapping={'file_name': 'file.name',
                       'covid': 'COVID_status',
                       'dataset': 'folder'},
        dataset='fmba_2021',
        get_extra_info=False,
        chain_to_read='TRB',
        functional_seqs_type='nonfunctional'
    ).to_csv('../data/usage_matrix_nonfunctional.csv')


if __name__ == "__main__":
    if 'snakemake' in globals():
        starting_pattern = 'TR' + ('B' if snakemake.params.gene == 'TRB' else 'A') + (
            'V' if snakemake.params.gene_type == 'v' else 'J')
        if snakemake.params.gene == 'TRB':
            if snakemake.params.platform == 'fmba':
                run_beta_chain_fmba_matrix_creation(desc_path=snakemake.params.desc_path,
                                                    output_path=snakemake.output[0],
                                                    raw_data_path=snakemake.params.raw_data_path,
                                                    gene_type=snakemake.params.gene_type)
            if snakemake.params.platform == 'adaptive':
                adaptive_hip_matrix_creation(hip_desc_path=snakemake.params.hip_desc_path,
                                             adaptive_desc_path=snakemake.params.adaptive_desc_path,
                                             raw_path=snakemake.params.raw_data_path,
                                             save_path=snakemake.output[0],
                                             gene_type=snakemake.params.gene_type)

        if snakemake.params.gene == 'TRA':
            if snakemake.params.platform == 'fmba':
                run_alpha_chain_fmba_matrix_creation(desc_path=snakemake.params.desc_path,
                                                     output_path=snakemake.output[0],
                                                     raw_data_path=snakemake.params.raw_data_path,
                                                     gene_type=snakemake.params.gene_type)
