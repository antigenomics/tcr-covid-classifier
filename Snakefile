configfile: 'config.yaml'

rule fmba_beta_usage_matrix_creation:
    threads: 1
    params: desc_path=f'{config["fmba_txt_desc"]}',
            raw_data_path=f'{config["all_raw_data_path"]}',
            gene='TRB',
            platform='fmba'
    output: 'data/usage_matrix_fmba_TRB.csv'
    script: 'source/repertoire_matrix_extraction.py'

rule fmba_alpha_usage_matrix_creation:
    threads: 1
    params: desc_path=f'{config["fmba_txt_desc"]}',
            raw_data_path=f'{config["all_raw_data_path"]}',
            gene='TRA',
            platform='fmba'
    output: 'data/usage_matrix_fmba_TRA.csv'
    script: 'source/repertoire_matrix_extraction.py'

rule fmba_beta_usage_matrix_standardization:
    threads: 1
    input: 'data/usage_matrix_fmba_TRB.csv'
    params: gene='TRB',
            platform='fmba'
    output: 'data/standardized_usage_matrix_fmba_TRB.csv', 'data/normalized_usage_matrix_fmba_TRB.csv'
    script: 'source/usage_matrix_preprocessing.py'

rule fmba_alpha_usage_matrix_standardization:
    threads: 1
    input: 'data/usage_matrix_fmba_TRA.csv'
    params: gene='TRA',
            platform='fmba'
    output: 'data/standardized_usage_matrix_fmba_TRA.csv', 'data/normalized_usage_matrix_fmba_TRA.csv'
    script: 'source/usage_matrix_preprocessing.py'

rule figure_1:
    threads: 1
    input: 'data/desc_fmba_not_nan_hla.csv', 'data/fmba_2021.txt', 'data/desc_fmba_new_split.csv',
         'data/standardized_usage_matrix_fmba_TRB.csv', 'data/standardized_usage_matrix_fmba_TRB.csv',
         'data/standardized_usage_matrix_fmba_TRA.csv', 'data/normalized_usage_matrix_fmba_TRA.csv'
    output: 'figures/fig1.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/fig1.ipynb
            python publication-notebooks/fig1.py
           '''

rule resampling_fmba_beta:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRB.csv'
    params: gene='TRB',
            platform='fmba'
    output: directory(f'{config["all_raw_data_path"]}/downsampled_fmba_TRB')
    script: 'source/downsampling_procedure.py'

rule resampling_fmba_alpha:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRA.csv'
    params: gene='TRA',
            platform='fmba'
    output: directory(f'{config["all_raw_data_path"]}/downsampled_fmba_TRA')
    script: 'source/downsampling_procedure.py'

rule clones_extraction_fmba_beta:
    threads: 2
    resources: mem="10GB"
    input: 'data/standardized_usage_matrix_fmba_TRB.csv', f'{config["all_raw_data_path"]}/downsampled_fmba_TRB'
    params: n_clones=500000,
            platform='fmba'
    output: 'data/most_used_500k_fmba_TRB.csv'
    script: 'source/clonotypes_extraction_procedure.py'

rule clones_extraction_fmba_alpha:
    threads: 2
    resources: mem="10GB"
    input: 'data/standardized_usage_matrix_fmba_TRA.csv', f'{config["all_raw_data_path"]}/downsampled_fmba_TRA'
    params: n_clones=500000,
            platform='fmba'
    output: 'data/most_used_500k_fmba_TRA.csv'
    script: 'source/clonotypes_extraction_procedure.py'

rule clone_matrix_creation_fmba_TRB:
    threads: 40
    resources: mem="10GB"
    params: platform='fmba'
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRB',
           'data/most_used_500k_fmba_TRB.csv'
    output: 'data/clone_matrix_fmba_TRB_top_500k.csv'
    script: 'source/clonotypes_matrix_creation.py'

rule create_sample_to_num_of_clones_mapping_fmba_TRB:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRB',
    output: 'data/run_to_number_of_clones_fmba_TRB.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

rule clone_matrix_creation_fmba_TRA:
    threads: 40
    resources: mem="10GB"
    params: platform='fmba'
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRA',
           'data/most_used_500k_fmba_TRA.csv'
    output: 'data/clone_matrix_fmba_TRA_top_500k.csv'
    script: 'source/clonotypes_matrix_creation.py'

rule create_sample_to_num_of_clones_mapping_fmba_TRA:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRA',
    output: 'data/run_to_number_of_clones_fmba_TRA.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

rule fisher_test_fmba_TRA:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           'data/run_to_number_of_clones_fmba_TRA.csv',
           'data/clone_matrix_fmba_TRA_top_500k.csv'
    params: n_clones=500000,
            platform='fmba'
    output: 'data/covid_significant_clones_fmba_TRA_top_500k.csv',
            'data/covid_significant_clone_pvals_fmba_TRA_top_500k.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule fisher_test_fmba_TRB:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           'data/run_to_number_of_clones_fmba_TRB.csv',
           'data/clone_matrix_fmba_TRB_top_500k.csv'
    params: n_clones=500000,
            platform='fmba'
    output: 'data/covid_significant_clones_fmba_TRB_top_500k.csv',
            'data/covid_significant_clone_pvals_fmba_TRB_top_500k.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule fisher_significant_clone_matrix_fmba_TRB:
    input: 'data/clone_matrix_fmba_TRB_top_500k.csv',
           'data/covid_significant_clones_fmba_TRB_top_500k.csv'
    output: 'data/significant_clone_matrix_fisher_fmba_TRB_top_500k.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule fisher_significant_clone_matrix_fmba_TRA:
    input: 'data/clone_matrix_fmba_TRA_top_500k.csv',
           'data/covid_significant_clones_fmba_TRA_top_500k.csv'
    output: 'data/significant_clone_matrix_fisher_fmba_TRA_top_500k.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule olga_pgen_generation_fmba_TRB:
    input: 'data/covid_significant_clones_fmba_TRB_top_500k.csv'
    output: 'data/covid_fmba_TRB_pgen.csv'
    shell: 'olga-compute_pgen --humanTRB -i {input} > {output}'

rule olga_pgen_generation_fmba_TRA:
    input: 'data/covid_significant_clones_fmba_TRA_top_500k.csv'
    output: 'data/covid_fmba_TRA_pgen.csv'
    shell: 'olga-compute_pgen --humanTRB -i {input} > {output}'