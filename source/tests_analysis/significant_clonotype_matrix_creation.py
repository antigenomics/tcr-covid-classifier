import pandas as pd


def create_significant_clonotype_matrix(clonotype_matrix_path, significant_clones_path, save_path):
    clonotype_matrix = pd.read_csv(clonotype_matrix_path)
    clonotype_matrix = clonotype_matrix.set_index('cdr3aa').T.reset_index().rename(columns={'index': 'run'})
    significant_clones = pd.read_csv(significant_clones_path)
    clone_column = 'clone'
    if clone_column not in significant_clones.columns:
        clone_column = 'cdr3aa'
    useful_columns = [x for x in significant_clones[clone_column] if x in clonotype_matrix.columns] + ['run']
    clonotype_matrix[useful_columns].to_csv(save_path)


def sign_clone_matrices_for_all_alleles_covid():
    hla_keys = pd.read_csv('data/hla_keys.csv')['0']
    print(hla_keys)
    for hla in list(hla_keys)[6:]:
        print(f'Started processing {hla}')
        significant_clones_path = f'data/hla_covid_results/covid_clones_500k_top_1_mismatch_hla_{hla}.csv'
        clonotype_matrix_path = f'data/hla_clonotype_matrix/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{hla}.csv'
        save_path = f'data/hla_sign_clone_matrix/hla_covid_clonotype_matrix_500k_top_1_mismatch_hla_{hla}.csv'
        create_significant_clonotype_matrix(clonotype_matrix_path=clonotype_matrix_path,
                                            significant_clones_path=significant_clones_path,
                                            save_path=save_path)


def sign_clone_matrices_hla():
    hla_keys = pd.read_csv('data/hla_keys.csv')['0']
    print(hla_keys)
    for hla in list(hla_keys):
        print(f'Started processing {hla}')
        try:
            significant_clones_path = f'data/hla_associated_clones_fisher/hla_associated_clones_500k_top_1_mismatch_hla_{hla}.csv'
            clonotype_matrix_path = f'data/hla_clonotype_matrix/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{hla}.csv'
            save_path = f'data/hla_associated_sign_clone_matrix/hla_clonotype_matrix_500k_top_1_mismatch_hla_{hla}.csv'
            create_significant_clonotype_matrix(clonotype_matrix_path=clonotype_matrix_path,
                                                significant_clones_path=significant_clones_path,
                                                save_path=save_path)
        except Exception as e:
            print('This allele is not ready yet')


def sign_clone_matrices_hla_with_adaptive():
    hla_keys = pd.read_csv('data/hla_keys.csv')['0']
    print(hla_keys)
    for hla in list(hla_keys):
        print(f'Started processing {hla}')
        # try:
        significant_clones_path = f'data/hla_associated_clones_fisher/hla_associated_clones_500k_top_1_mismatch_hla_{hla}.csv'
        clonotype_matrix_path = f'data/hla_sign_clone_matrix/clones_for_all_alleles_matrix_500k_top_1_mismatch.csv'
        save_path = f'data/hla_associated_sign_clone_matrix/hla_clonotype_matrix_500k_top_1_mismatch_hla_{hla}_with_adaptive.csv'
        create_significant_clonotype_matrix(clonotype_matrix_path=clonotype_matrix_path,
                                            significant_clones_path=significant_clones_path,
                                            save_path=save_path)
        # except Exception as e:
        #     print('This allele is not ready yet')


def sign_clone_matrices_for_all_alleles_covid_based():
    hla_keys = pd.read_csv('data/hla_keys.csv')['0']
    print(hla_keys)
    for hla in list(hla_keys):
        print(f'Started processing {hla}')
        clonotype_matrix_path = f'data/covid_clonotype_matrix/clonotype_matrix_covid_fmba_top_500k_1_mismatch.csv'
        significant_clones_path = f'data/covid_associated_clones/significant_clones_500k_top_1_mismatch_hla_{hla}.csv'
        save_path = f'data/covid_associated_clone_matrix/hla_covid_clonotype_matrix_500k_top_1_mismatch_hla_{hla}.csv'
        try:
            create_significant_clonotype_matrix(clonotype_matrix_path=clonotype_matrix_path,
                                                significant_clones_path=significant_clones_path,
                                                save_path=save_path)
        except Exception:
            print(f'No info for allele {hla}')


if __name__ == "__main__":
    if 'snakemake' in globals():
        if snakemake.params.platform == 'fmba':
            create_significant_clonotype_matrix(clonotype_matrix_path=snakemake.input[0],
                                                significant_clones_path=snakemake.input[1],
                                                save_path=snakemake.output[0])
        if snakemake.params.platform == 'allele':
            if snakemake.params.hla_to_consider == []:
                hla_keys = snakemake.params.hla_to_consider
            else:
                hla_keys = pd.read_csv('data/hla_keys.csv')['0']
            print(hla_keys)
            import os

            if not os.path.exists(snakemake.output[0]):
                os.mkdir(snakemake.output[0])
            for hla in hla_keys:
                create_significant_clonotype_matrix(
                    clonotype_matrix_path=f'{snakemake.input[1]}/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{hla}.csv',
                    significant_clones_path=f'{snakemake.input[0]}/hla_associated_clones_500k_top_1_mismatch_hla_{hla}.csv',
                    save_path=f'{snakemake.output[0]}/hla_covid_clonotype_matrix_500k_top_1_mismatch_hla_{hla}.csv'
                )
        if snakemake.params.platform == 'fmba-allele':
            if snakemake.params.hla_to_consider == []:
                hla_keys = snakemake.params.hla_to_consider
            else:
                hla_keys = pd.read_csv('data/hla_keys.csv')['0']
            print(hla_keys)
            import os
            if not os.path.exists(snakemake.output[0]):
                os.mkdir(snakemake.output[0])
            for hla in hla_keys:
                joint_clones = pd.read_csv(
                    f'{snakemake.input[0]}/hla_associated_clones_500k_top_1_mismatch_hla_{hla}.csv').merge(
                    pd.read_csv(
                        f'{snakemake.input[1]}/hla_covid_associated_clones_500k_top_1_mismatch_hla_{hla}.csv'
                    )
                ).to_csv(
                    f'{snakemake.output[0]}/hla_covid_joint_sign_clones_{hla}.csv', index=False
                )

                create_significant_clonotype_matrix(
                    clonotype_matrix_path=f'{snakemake.input[2]}/clonotype_matrix_500k_1_mismatch_top_fmba_hla_{hla}.csv',
                    significant_clones_path=f'{snakemake.output[0]}/hla_covid_joint_sign_clones_{hla}.csv',
                    save_path=f'{snakemake.output[0]}/hla_covid_clonotype_matrix_500k_top_1_mismatch_hla_{hla}.csv'
                )
