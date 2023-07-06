configfile: 'config.yaml'

rule all:
  input:
    'figures/fig1.png',
    'figures/fig2.png',
    'figures/fig3.png',
    'figures/fig4.png',
    'figures/fig5.png',
    'figures/supp_fig2.png',
    'figures/supp_fig3.png'

rule preprocess_metadata:
    input: 'data/fmba_2021.txt'
    output: 'data/preprocessed_fmba_metadata.csv'
    script: 'source/preprocess_raw_desc_file.py'

rule fmba_beta_usage_matrix_creation:
    threads: 1
    input: f'{config["fmba_desc"]}'
    params: desc_path=f'{config["fmba_desc"]}',
            raw_data_path=f'{config["all_raw_data_path"]}',
            gene='TRB',
            platform='fmba'
    output: 'data/usage_matrix_fmba_TRB.csv'
    script: 'source/repertoire_matrix_extraction.py'

rule fmba_beta_usage_matrix_standardization:
    threads: 1
    input: 'data/usage_matrix_fmba_TRB.csv'
    params: gene='TRB',
            platform='fmba'
    output: 'data/standardized_usage_matrix_fmba_TRB_wo_test_runs.csv', 'data/normalized_usage_matrix_fmba_TRB.csv'
    script: 'source/usage_matrix_preprocessing.py'

rule fmba_beta_usage_matrix_test_column_creation:
    threads: 1
    input: 'data/standardized_usage_matrix_fmba_TRB_wo_test_runs.csv'
    params: test_selection_method='batch', test_batch='NovaSeq6'
    output: 'data/standardized_usage_matrix_fmba_TRB.csv'
    script: 'source/select_test_runs.py'

rule resampling_fmba_beta:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRB.csv'
    params: gene='TRB',
            platform='fmba'
    output: directory(f'{config["all_raw_data_path"]}/downsampled_fmba_TRB')
    script: 'source/downsampling_procedure.py'

rule clones_extraction_fmba_beta:
    threads: 2
    resources: mem="10GB"
    input: 'data/standardized_usage_matrix_fmba_TRB.csv', f'{config["all_raw_data_path"]}/downsampled_fmba_TRB'
    params: n_clones=2,
            platform='fmba',
            sampling_method='unique-occurence'
    output: 'data/most_used_500k_fmba_TRB.csv'
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

rule clone_matrix_creation_vdjdb_fmba_TRB:
    threads: 40
    resources: mem="10GB"
    params: platform='vdjdb'
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRB',
           'data/vdjdb_TRB_covid_clones.csv'
    output: 'data/clone_matrix_vdjdb_fmba_TRB_not_processed.csv'
    script: 'source/clonotypes_matrix_creation.py'

rule fisher_significant_clone_matrix_vdjdb_fmba_TRB:
    input: 'data/clone_matrix_vdjdb_fmba_TRB_not_processed.csv',
           'data/vdjdb_TRB_covid_clones.csv'
    params: platform='fmba'
    output: 'data/significant_clone_matrix_fisher_fmba_TRB_vdjdb.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule create_sample_to_num_of_clones_mapping_fmba_TRB:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRB',
    params: make_read_count_col=False
    output: 'data/run_to_number_of_clones_fmba_TRB.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

rule fisher_test_fmba_TRB:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           'data/run_to_number_of_clones_fmba_TRB.csv',
           'data/clone_matrix_fmba_TRB_top_500k.csv'
    params: platform='fmba', significant_threshold=0.01, drop_test=False
    output: 'data/covid_significant_clones_fmba_TRB_top_500k.csv',
            'data/covid_significant_clone_pvals_fmba_TRB_top_500k.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule fisher_test_fmba_TRB_random:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           'data/run_to_number_of_clones_fmba_TRB.csv',
           'data/clone_matrix_fmba_TRB_top_500k.csv'
    params: platform='fmba-random-488', significant_threshold=0.01, drop_test=False
    output: 'data/covid_significant_clones_fmba_TRB_random_488.csv',
            'data/covid_significant_clone_pvals_fmba_TRB_random_488.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule fisher_significant_clone_matrix_fmba_TRB_random:
    input: 'data/clone_matrix_fmba_TRB_top_500k.csv',
           'data/covid_significant_clones_fmba_TRB_random_488.csv'
    params: platform='fmba'
    output: 'data/significant_clone_matrix_fisher_fmba_TRB_random_488.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule fisher_test_fmba_TRB_vdjdb:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           'data/run_to_number_of_clones_fmba_TRB.csv',
           'data/clone_matrix_vdjdb_fmba_TRB_not_processed.csv'
    params: platform='vdjdb', significant_threshold=1, drop_test=True
    output: 'data/covid_significant_clones_fmba_TRB_vdjdb.csv',
            'data/covid_significant_clone_pvals_fmba_TRB_vdjdb.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule fisher_significant_clone_matrix_fmba_TRB:
    input:  'data/clone_matrix_fmba_TRB_top_500k.csv',
           'data/covid_significant_clones_fmba_TRB_top_500k.csv'
    params: platform='fmba'
    output: 'data/significant_clone_matrix_fisher_fmba_TRB_top_500k.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule olga_pgen_generation_fmba_TRB:
    input: 'data/covid_significant_clones_fmba_TRB_top_500k.csv'
    output: 'data/covid_fmba_TRB_pgen.csv'
    shell: 'olga-compute_pgen --humanTRB -i {input} > {output}'

rule fisher_significant_clone_matrix_wo_leaks_fmba_TRB:
    input: 'data/clone_matrix_fmba_TRB_top_500k.csv',
            'data/covid_fmba_TRB_pgen.csv',
            'data/run_to_number_of_clones_fmba_TRB.csv'
    params: pgen_threshold=1e-9, preal_threshold=1e-6
    output: 'data/significant_clone_matrix_fisher_fmba_TRB_top_500k_wo_leaks.csv',
    script: 'source/leaks_deletion.py'

rule create_cluster_bool_matrix_TRB:
    input: 'data/significant_clone_matrix_fisher_fmba_TRB_top_500k_wo_leaks.csv'
    params: gene='TRB'
    output: 'data/cluster_presence_matrix_TRB.csv'
    script: 'source/make_cluster_usage_table.py'
########################################################################################################################
rule create_hla_desc_files:
    input: 'data/desc_fmba_not_nan_hla.csv'
    output: 'data/hla_keys.csv', directory('data/hla_desc')
    script: 'source/create_hla_description_files.py'

rule hla_specific_clones_extraction_TRB:
    threads: 2
    resources: mem="10GB"
    input: 'data/hla_desc', 'data/standardized_usage_matrix_fmba_TRB.csv', f'{config["all_raw_data_path"]}/downsampled_fmba_TRB'
    params: platform='fmba-allele',
            hla_to_consider=['A*02', 'DQB1*05', 'DRB1*16']
    output: directory('data/hla_most_used_clones')
    script: 'source/clonotypes_extraction_procedure.py'


rule hla_specific_clones_extraction_TRA:
    threads: 2
    resources: mem="10GB"
    input: 'data/hla_desc', 'data/standardized_usage_matrix_fmba_TRA.csv', f'{config["all_raw_data_path"]}/downsampled_fmba_TRA'
    params: platform='fmba-allele',
            hla_to_consider=['A*02']
    output: directory('data/hla_most_used_clones_TRA')
    script: 'source/clonotypes_extraction_procedure.py'

rule hla_clone_matrix_creation_TRB:
    threads: 40
    resources: mem="10GB"
    params: platform='fmba-allele',
            hla_to_consider=['A*02', 'DQB1*05', 'DRB1*16']
    input: 'data/standardized_usage_matrix_fmba_TRB.csv',
           'data/hla_most_used_clones',
            f'{config["all_raw_data_path"]}/downsampled_fmba_TRB'
    output: directory('data/hla_clonotype_matrices')
    script: 'source/clonotypes_matrix_creation.py'

rule hla_clone_matrix_creation_TRA:
    threads: 40
    resources: mem="10GB"
    params: platform='fmba-allele',
            hla_to_consider=['A*02']
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           'data/hla_most_used_clones_TRA',
            f'{config["all_raw_data_path"]}/downsampled_fmba_TRA'
    output: directory('data/hla_clonotype_matrices_TRA')
    script: 'source/clonotypes_matrix_creation.py'

rule hla_associatiated_TRB_clones_search:
    threads: 40
    params: hla_to_consider=['A*02', 'DQB1*05', 'DRB1*16']
    input: 'data/run_to_number_of_clones_fmba_TRB.csv',
            'data/hla_desc',
            'data/hla_clonotype_matrices',
    output: directory('data/hla_associated_clones')
    script: 'source/tests_analysis/hla_test_new.py'

rule hla_based_covid_associated_TRB_clones_search:
    threads: 40
    params: hla_to_consider=['A*02', 'DQB1*05', 'DRB1*16'],
            platform='fmba-allele', significant_threshold=0.01, drop_test=False
    input: 'data/run_to_number_of_clones_fmba_TRB.csv',
            'data/hla_desc',
            'data/hla_clonotype_matrices',
    output: directory('data/hla_covid_associated_clones')
    script: 'source/tests_analysis/covid_test_250k.py'

rule hla_based_covid_associated_TRA_clones_search:
    threads: 40
    params: hla_to_consider=['A*02'],
            platform='fmba-allele', significant_threshold=0.01, drop_test=False
    input: 'data/run_to_number_of_clones_fmba_TRA.csv',
            'data/hla_desc',
            'data/hla_clonotype_matrices_TRA',
    output: directory('data/hla_covid_associated_clones_TRA')
    script: 'source/tests_analysis/covid_test_250k.py'

rule hla_covid_significant_usage_matrix_creation_TRB:
    input:  'data/hla_associated_clones',
            'data/hla_covid_associated_clones',
            'data/hla_clonotype_matrices',
    params: platform='fmba-allele',
            hla_to_consider=['A*02', 'DQB1*05', 'DRB1*16']
    output: directory('data/hla_sign_clone_matrix')
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule create_a02_significant_covid_clone_matrix_TRB:
    input:  'data/clone_matrix_fmba_TRB_top_500k.csv',
           'data/hla_covid_associated_clones/hla_covid_associated_clones_500k_top_1_mismatch_hla_A*02.csv'
    params: platform='fmba'
    output: 'data/significant_clone_matrix_fisher_fmba_TRB_A02.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule create_a02_significant_covid_clone_matrix_TRA:
    input:  'data/clone_matrix_fmba_TRA_top_500k.csv',
           'data/hla_covid_associated_clones_TRA/hla_covid_associated_clones_500k_top_1_mismatch_hla_A*02.csv'
    params: platform='fmba'
    output: 'data/significant_clone_matrix_fisher_fmba_TRA_A02.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

########################################################################################################################

rule fmba_alpha_usage_matrix_standardization:
    threads: 1
    input: 'data/usage_matrix_fmba_TRA.csv'
    params: gene='TRA',
            platform='fmba'
    output: 'data/standardized_usage_matrix_fmba_TRA_wo_test_runs.csv', 'data/normalized_usage_matrix_fmba_TRA.csv'
    script: 'source/usage_matrix_preprocessing.py'

rule fmba_alpha_usage_matrix_test_column_creation:
    threads: 1
    input: 'data/standardized_usage_matrix_fmba_TRA_wo_test_runs.csv'
    params: test_selection_method='batch', test_batch='NovaSeq6'
    output: 'data/standardized_usage_matrix_fmba_TRA.csv'
    script: 'source/select_test_runs.py'

rule fmba_alpha_usage_matrix_creation:
    threads: 1
    input: f'{config["fmba_desc"]}'
    params: desc_path=f'{config["fmba_desc"]}',
            raw_data_path=f'{config["all_raw_data_path"]}',
            gene='TRA',
            platform='fmba'
    output: 'data/usage_matrix_fmba_TRA.csv'
    script: 'source/repertoire_matrix_extraction.py'

rule resampling_fmba_alpha:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRA.csv'
    params: gene='TRA',
            platform='fmba'
    output: directory(f'{config["all_raw_data_path"]}/downsampled_fmba_TRA')
    script: 'source/downsampling_procedure.py'

rule clones_extraction_fmba_alpha:
    threads: 2
    resources: mem="10GB"
    input: 'data/standardized_usage_matrix_fmba_TRA.csv', f'/projects/fmba_covid/1_data_links/downsampled_alpha'
    params: n_clones=2,
            platform='fmba',
            sampling_method='unique-occurence'
    output: 'data/most_used_500k_fmba_TRA.csv'
    script: 'source/clonotypes_extraction_procedure.py'

rule clone_matrix_creation_fmba_TRA:
    threads: 40
    resources: mem="10GB"
    params: platform='fmba'
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRA',
           'data/most_used_500k_fmba_TRA.csv'
    output: 'data/clone_matrix_fmba_TRA_top_500k.csv'
    script: 'source/clonotypes_matrix_creation.py'

rule clone_matrix_creation_vdjdb_fmba_TRA:
    threads: 40
    resources: mem="10GB"
    params: platform='vdjdb'
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRA',
           'data/vdjdb_TRA_covid_clones.csv'
    output: 'data/clone_matrix_vdjdb_fmba_TRA_not_processed.csv'
    script: 'source/clonotypes_matrix_creation.py'

rule fisher_significant_clone_matrix_vdjdb_fmba_TRA:
    input: 'data/clone_matrix_vdjdb_fmba_TRA_not_processed.csv',
           'data/vdjdb_TRA_covid_clones.csv'
    params: platform='fmba'
    output: 'data/significant_clone_matrix_fisher_fmba_TRA_vdjdb.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule create_sample_to_num_of_clones_mapping_fmba_TRA:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           f'{config["all_raw_data_path"]}/downsampled_fmba_TRA',
    params: make_read_count_col=False
    output: 'data/run_to_number_of_clones_fmba_TRA.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

rule fisher_test_fmba_TRA:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           'data/run_to_number_of_clones_fmba_TRA.csv',
           'data/clone_matrix_fmba_TRA_top_500k.csv'
    params: platform='fmba', significant_threshold=0.01, drop_test=False
    output: 'data/covid_significant_clones_fmba_TRA_top_500k.csv',
            'data/covid_significant_clone_pvals_fmba_TRA_top_500k.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule fisher_test_fmba_TRA_random:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           'data/run_to_number_of_clones_fmba_TRA.csv',
           'data/clone_matrix_fmba_TRA_top_500k.csv'
    params: platform='fmba-random-488', significant_threshold=0.01, drop_test=False
    output: 'data/covid_significant_clones_fmba_TRA_random_488.csv',
            'data/covid_significant_clone_pvals_fmba_TRA_random_488.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule fisher_test_fmba_TRA_vdjdb:
    threads: 40
    input: 'data/standardized_usage_matrix_fmba_TRA.csv',
           'data/run_to_number_of_clones_fmba_TRA.csv',
           'data/clone_matrix_vdjdb_fmba_TRA_not_processed.csv'
    params: platform='vdjdb', significant_threshold=1, drop_test=True
    output: 'data/covid_significant_clones_fmba_TRA_vdjdb.csv',
            'data/covid_significant_clone_pvals_fmba_TRA_vdjdb.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule fisher_significant_clone_matrix_fmba_TRA:
    input: 'data/clone_matrix_fmba_TRA_top_500k.csv',
           'data/covid_significant_clones_fmba_TRA_top_500k.csv'
    params: platform='fmba'
    output: 'data/significant_clone_matrix_fisher_fmba_TRA_top_500k.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule fisher_significant_clone_matrix_fmba_TRA_random:
    input: 'data/clone_matrix_fmba_TRA_top_500k.csv',
           'data/covid_significant_clones_fmba_TRA_random_488.csv'
    params: platform='fmba'
    output: 'data/significant_clone_matrix_fisher_fmba_TRA_random_488.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule olga_pgen_generation_fmba_TRA:
    input: 'data/covid_significant_clones_fmba_TRA_top_500k.csv'
    output: 'data/covid_fmba_TRA_pgen.csv'
    shell: 'olga-compute_pgen --humanTRA -i {input} > {output}'

rule fisher_significant_clone_matrix_wo_leaks_fmba_TRA:
    input: 'data/significant_clone_matrix_fisher_fmba_TRA_top_500k.csv',
            'data/covid_fmba_TRA_pgen.csv',
            'data/run_to_number_of_clones_fmba_TRA.csv'
    params: pgen_threshold=1e-9, preal_threshold=1e-6
    output: 'data/significant_clone_matrix_fisher_fmba_TRA_top_500k_wo_leaks.csv',
    script: 'source/leaks_deletion.py'

rule create_cluster_bool_matrix_TRA:
    input: 'data/significant_clone_matrix_fisher_fmba_TRA_top_500k_wo_leaks.csv'
    params: gene='TRA'
    output: 'data/cluster_presence_matrix_TRA.csv'
    script: 'source/make_cluster_usage_table.py'

rule TRA_TRB_pairing_analysis:
    input: 'data/significant_clone_matrix_fisher_fmba_TRA_top_500k_wo_leaks.csv',
           'data/significant_clone_matrix_fisher_fmba_TRB_top_500k_wo_leaks.csv',
            'data/vdjdb.txt',
            'data/cluster_presence_matrix_TRB.csv',
            'data/cluster_presence_matrix_TRA.csv'
    output: 'data/TRA_TRB_cooccurence_matrix_cooccurence_85.csv', 'data/alpha_beta_paired_epitopes.csv',
            'data/clone_matrix_covid_fmba_TRA_metaclone.csv', 'data/clone_matrix_covid_fmba_TRB_metaclone.csv',
            'figures/cooccured_epitopes_fmba.pdf'
    script: 'source/alpha_beta_paired_clones_search.py'

######################################################################################################################
rule adaptive_usage_matrix_creation:
    threads: 1
    params: hip_desc_path=f'{config["hip_desc"]}',
            adaptive_desc_path=f'{config["adaptive_desc"]}',
            raw_data_path=f'{config["all_raw_data_path"]}',
            gene='TRB',
            platform='adaptive'
    output: 'data/usage_matrix_adaptive.csv'
    script: 'source/repertoire_matrix_extraction.py'

rule adaptive_usage_matrix_standardization:
    threads: 1
    input: 'data/usage_matrix_adaptive.csv'
    params: gene='TRB',
            platform='adaptive'
    output: 'data/standardized_usage_matrix_adaptive.csv', 'data/normalized_usage_matrix_adaptive.csv'
    script: 'source/usage_matrix_preprocessing.py'

rule resampling_adaptive:
    threads: 40
    input: 'data/standardized_usage_matrix_adaptive.csv'
    params: gene='TRB',
            platform='adaptive'
    output: directory(f'{config["all_raw_data_path"]}/downsampled_adaptive_TRB')
    script: 'source/downsampling_procedure.py'

rule clones_extraction_adaptive:
    threads: 2
    resources: mem="50GB"
    input: 'data/standardized_usage_matrix_adaptive.csv', f'{config["all_raw_data_path"]}/downsampled_adaptive_TRB'
    params: n_clones=500000,
            platform='adaptive',
            sampling_method='top'
    output: 'data/most_used_500k_adaptive.csv'
    script: 'source/clonotypes_extraction_procedure.py'

rule clone_matrix_creation_adaptive:
    threads: 40
    params: platform='adaptive'
    input: 'data/standardized_usage_matrix_adaptive.csv',
           f'{config["all_raw_data_path"]}/downsampled_adaptive_TRB',
           'data/most_used_500k_adaptive.csv'
    output: 'data/clone_matrix_adaptive_top_500k.csv'
    script: 'source/clonotypes_matrix_creation.py'

rule create_sample_to_num_of_clones_mapping_adaptive:
    threads: 40
    input: 'data/standardized_usage_matrix_adaptive.csv',
           f'{config["all_raw_data_path"]}/downsampled_adaptive_TRB',
    params: make_read_count_col=False
    output: 'data/run_to_number_of_clones_adaptive.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

rule fisher_test_adaptive:
    threads: 20
    input: 'data/standardized_usage_matrix_adaptive.csv',
           'data/run_to_number_of_clones_adaptive.csv',
           'data/clone_matrix_adaptive_top_500k.csv'
    params: platform='adaptive', significant_threshold=0.01 / 1e6, drop_test=False
    output: 'data/covid_significant_clones_adaptive.csv',
            'data/covid_significant_clone_pvals_adaptive.csv'
    script: 'source/tests_analysis/covid_test_250k.py'

rule joint_usage_matrix_creation:
    threads: 1
    input: 'data/usage_matrix_adaptive.csv',
            'data/usage_matrix_fmba_TRB.csv'
    params: gene='TRB',
            platform='joint'
    output: 'data/standardized_usage_matrix_joint.csv', 'data/normalized_usage_matrix_joint.csv'
    script: 'source/usage_matrix_preprocessing.py'

rule joint_clonotypes_resampling:
    threads: 40
    input: 'data/standardized_usage_matrix_joint.csv'
    params: gene='TRB',
            platform='joint'
    output: directory(f'{config["all_raw_data_path"]}/downsampled_joint')
    script: 'source/downsampling_procedure.py'

rule joint_clonotypes_resampling_50k:
    threads: 40
    input: 'data/standardized_usage_matrix_joint.csv'
    params: gene='TRB',
            platform='joint_50k'
    output: directory(f'{config["all_raw_data_path"]}/downsampled_joint_50k')
    script: 'source/downsampling_procedure.py'

rule clone_matrix_creation_joint_fmba_based:
    threads: 40
    resources: mem="10GB"
    params: platform='fmba'
    input: 'data/standardized_usage_matrix_joint.csv',
           f'{config["all_raw_data_path"]}/downsampled_joint',
           'data/covid_significant_clones_fmba_TRB_top_500k.csv'
    output: 'data/clone_matrix_joint_fmba_based.csv'
    script: 'source/clonotypes_matrix_creation.py'

rule clone_matrix_prep_joint_fmba_based:
    input:  'data/clone_matrix_joint_fmba_based.csv',
           'data/covid_significant_clones_fmba_TRB_top_500k.csv'
    params: platform='fmba'
    output: 'data/sign_clone_matrix_joint_fmba_based.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule clone_matrix_creation_joint_adaptive_based:
    threads: 40
    resources: mem="10GB"
    params: platform='adaptive'
    input: 'data/standardized_usage_matrix_joint.csv',
           f'{config["all_raw_data_path"]}/downsampled_joint_50k',
           'data/covid_significant_clones_adaptive.csv'
    output: 'data/clone_matrix_joint_adaptive_based.csv'
    script: 'source/clonotypes_matrix_creation.py'

rule clone_matrix_prep_joint_adaptive_based:
    input:  'data/clone_matrix_joint_adaptive_based.csv',
           'data/covid_significant_clones_adaptive.csv'
    params: platform='fmba'
    output: 'data/sign_clone_matrix_joint_adaptive_based.csv'
    script: 'source/tests_analysis/significant_clonotype_matrix_creation.py'

rule olga_pgen_generation_adaptive:
    input: 'data/covid_significant_clones_adaptive.csv'
    output: 'data/covid_adaptive_pgen.csv'
    shell: 'olga-compute_pgen --humanTRB -i {input} > {output}'

rule create_sample_to_num_of_clones_mapping_joint:
    threads: 40
    input: 'data/standardized_usage_matrix_joint.csv',
           f'{config["all_raw_data_path"]}/downsampled_joint',
    params: make_read_count_col=False
    output: 'data/run_to_number_of_clones_joint.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

rule create_sample_to_num_of_clones_mapping_joint_50k:
    threads: 40
    input: 'data/standardized_usage_matrix_joint.csv',
           f'{config["all_raw_data_path"]}/downsampled_joint_50k',
    params: make_read_count_col=False
    output: 'data/run_to_number_of_clones_joint_50k.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

######################################################################################################################

rule create_raw_data_reads_distribution_hip:
    threads: 40
    input: 'data/hip_full_prep.txt',
            f'{config["hip_raw_path"]}',
    params: make_read_count_col=True
    output: 'data/hip_with_reads_count.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

rule create_raw_data_reads_distribution_adaptive:
    threads: 40
    input: 'data/usage_matrix_adaptive.csv',
            f'{config["adaptive_raw_path"]}',
    params: make_read_count_col=True
    output: 'data/adaptive_with_reads_count.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

rule create_raw_data_reads_distribution_fmba:
    threads: 40
    input: '/projects/fmba_covid/fmba_upd_march/raw_filename_desc.csv',
           '/projects/fmba_covid/fmba_upd_march/data',
    params: make_read_count_col=True
    output: 'data/fmba_with_read_count.csv'
    script: 'source/get_sum_clonotypes_per_file.py'

######################################################################################################################

rule figure_1:
    threads: 1
    input: 'data/desc_fmba_not_nan_hla.csv', 'data/fmba_2021.txt', 'data/desc_fmba_new_split.csv',
         'data/standardized_usage_matrix_fmba_TRB.csv', 'data/normalized_usage_matrix_fmba_TRB.csv',
         'data/standardized_usage_matrix_fmba_TRA.csv', 'data/normalized_usage_matrix_fmba_TRA.csv',
         'publication-notebooks/fig1.ipynb'
    output: 'figures/fig1.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/fig1.ipynb
            python publication-notebooks/fig1.py
            rm publication-notebooks/fig1.py
           '''

rule figure_2:
    threads: 1
    input: 'data/vdjdb.txt', 'data/desc_fmba_not_nan_hla.csv',
           'data/significant_clone_matrix_fisher_fmba_TRB_top_500k.csv',
            'data/significant_clone_matrix_fisher_fmba_TRB_top_500k_wo_leaks.csv',
            'data/covid_fmba_TRB_pgen.csv',
            'data/run_to_number_of_clones_fmba_TRB.csv',
            'data/significant_clone_matrix_fisher_fmba_TRA_top_500k.csv',
            'data/significant_clone_matrix_fisher_fmba_TRA_top_500k_wo_leaks.csv',
            'data/covid_fmba_TRA_pgen.csv',
            'data/run_to_number_of_clones_fmba_TRA.csv',
            'data/TRA_TRB_cooccurence_matrix_cooccurence_85.csv',
            'data/alpha_beta_paired_epitopes.csv',
            'data/standardized_usage_matrix_fmba_TRB.csv',
            'data/standardized_usage_matrix_fmba_TRA.csv',
            'publication-notebooks/fig2.ipynb'
    output: 'figures/fig2.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/fig2.ipynb
            python publication-notebooks/fig2.py
            rm publication-notebooks/fig2.py
           '''

rule figure_3:
    threads: 1
    input: 'data/hla_sign_clone_matrix',
            'data/desc_fmba_not_nan_hla.csv',
            'data/run_to_number_of_clones_fmba_TRB.csv',
            'publication-notebooks/fig3.ipynb'
    output: 'figures/fig3.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/fig3.ipynb
            python publication-notebooks/fig3.py
            rm publication-notebooks/fig3.py
           '''

rule figure_4:
    threads: 1
    input: 'data/run_to_number_of_clones_fmba_TRA.csv',
            'data/run_to_number_of_clones_fmba_TRB.csv',
            'data/desc_fmba_not_nan_hla.csv',
            'data/significant_clone_matrix_fisher_fmba_TRB_top_500k_wo_leaks.csv',
            'data/significant_clone_matrix_fisher_fmba_TRA_top_500k_wo_leaks.csv',
            'data/hla_keys.csv',
            'data/clone_matrix_covid_fmba_TRB_metaclone.csv',
            'data/clone_matrix_covid_fmba_TRA_metaclone.csv',
            'data/hla_desc',
            'publication-notebooks/fig4.ipynb'
    output: 'figures/fig4.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/fig4.ipynb
            python publication-notebooks/fig4.py
            rm publication-notebooks/fig4.py
           '''

rule figure_5:
    threads: 1
    input: 'data/standardized_usage_matrix_joint.csv', 'data/normalized_usage_matrix_joint.csv',
         'data/covid_significant_clone_pvals_adaptive.csv', 'data/covid_adaptive_pgen.csv',
         'data/vdjdb.txt', 'data/run_to_number_of_clones_adaptive.csv', 'data/desc_fmba_not_nan_hla.csv',
         'data/run_to_number_of_clones_joint.csv', 'data/run_to_number_of_clones_fmba_TRB.csv',
         'data/sign_clone_matrix_joint_fmba_based.csv', 'data/sign_clone_matrix_joint_adaptive_based.csv',
         'data/run_to_number_of_clones_joint_50k.csv'
    output: 'figures/fig5.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/fig5.ipynb
            python publication-notebooks/fig5.py
            rm publication-notebooks/fig5.py
           '''

rule figure_6:
    threads: 1
    input: 'data/run_to_number_of_clones_fmba_TRA.csv',
            'data/standardized_usage_matrix_fmba_TRA.csv',
            'data/significant_clone_matrix_fisher_fmba_TRA_random_488.csv',
            'data/hla_keys.csv',
            'data/run_to_number_of_clones_fmba_TRB.csv',
            'data/significant_clone_matrix_fisher_fmba_TRB_random_488.csv',
            'data/hla_desc/fmba_desc_hla_A*02.csv',
            'data/significant_clone_matrix_fisher_fmba_TRB_A02.csv',
            'data/significant_clone_matrix_fisher_fmba_TRA_A02.csv',
    output: 'figures/fig6.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/fig6.ipynb
            python publication-notebooks/fig6.py
            rm publication-notebooks/fig6.py
           '''

rule supplementary_figure_2:
    threads: 1
    input: 'data/run_to_number_of_clones_fmba_TRA.csv',
            'data/standardized_usage_matrix_fmba_TRA.csv',
            'data/significant_clone_matrix_fisher_fmba_TRA_top_500k_wo_leaks.csv',
            'data/run_to_number_of_clones_fmba_TRB.csv',
            'data/hla_keys.csv',
            'data/clone_matrix_covid_fmba_TRB_metaclone.csv',
            'data/clone_matrix_covid_fmba_TRA_metaclone.csv',
    output: 'figures/supp_fig2.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/supp_fig2.ipynb
            python publication-notebooks/supp_fig2.py
            rm publication-notebooks/supp_fig2.py
           '''

rule supplementary_figure_3:
    threads: 1
    input: 'data/run_to_number_of_clones_fmba_TRA.csv',
            'data/run_to_number_of_clones_fmba_TRB.csv',
            'data/desc_fmba_not_nan_hla.csv',
            'data/significant_clone_matrix_fisher_fmba_TRB_vdjdb.csv',
            'data/significant_clone_matrix_fisher_fmba_TRA_vdjdb.csv',
            'data/hla_keys.csv', 'data/hla_desc',
            'publication-notebooks/supp_fig3.ipynb'
    output: 'figures/supp_fig3.png'
    shell: '''
            jupyter nbconvert --to python publication-notebooks/supp_fig3.ipynb
            python publication-notebooks/supp_fig3.py
            rm publication-notebooks/supp_fig3.py
           '''