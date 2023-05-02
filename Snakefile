configfile: 'config.yaml'

rule fmba_beta_usage_matrix_creation:
    params: desc_path=f'{config["fmba_txt_desc"]}',
            raw_data_path=f'{config["all_raw_data_path"]}',
            gene='TRB',
            platform='fmba'
    output: 'data/usage_matrix_fmba_TRB.csv'
    script: 'source/repertoire_matrix_extraction.py'

rule fmba_alpha_usage_matrix_creation:
    params: desc_path=f'{config["fmba_txt_desc"]}',
            raw_data_path=f'{config["all_raw_data_path"]}',
            gene='TRA',
            platform='fmba'
    output: 'data/usage_matrix_fmba_TRA.csv'
    script: 'source/repertoire_matrix_extraction.py'

rule fmba_beta_usage_matrix_standardization:
    input: 'data/usage_matrix_fmba_TRB.csv'
    params: gene='TRB',
            platform='fmba'
    output: 'data/standardized_usage_matrix_fmba_TRB.csv', 'data/normalized_usage_matrix_fmba_TRB.csv'
    script: 'source/usage_matrix_preprocessing.py'


rule fmba_alpha_usage_matrix_standardization:
    input: 'data/usage_matrix_fmba_TRA.csv'
    params: gene='TRA',
            platform='fmba'
    output: 'data/standardized_usage_matrix_fmba_TRA.csv', 'data/normalized_usage_matrix_fmba_TRA.csv'
    script: 'source/usage_matrix_preprocessing.py'

rule figure_1:
    input: 'data/desc_fmba_not_nan_hla.csv', 'data/fmba_2021.txt', 'data/desc_fmba_new_split.csv',
           'data/standardized_usage_matrix_fmba_TRB.csv', 'data/standardized_usage_matrix_fmba_TRB.csv',
           'data/standardized_usage_matrix_fmba_TRA.csv', 'data/normalized_usage_matrix_fmba_TRA.csv'
    output: 'figures/fig1.png'
    # conda: 'environment.yml'
    # notebook: 'publication-notebooks/fig1.ipynb'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/fig1.ipynb
            python publication-notebooks/fig1.py
           '''