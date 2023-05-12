import pandas as pd
from utils.data_utils import calculate_real_and_gen_proba, prepare_clonotype_matrix

if __name__ == "__main__":
    cm_path = snakemake.input[0]
    clone_pgen_path = snakemake.input[1]
    run_to_number_of_clones_path = snakemake.input[2]

    pgen_threshold = snakemake.params.pgen_threshold
    preal_threshold = snakemake.params.preal_threshold

    save_path = snakemake.output[0]

    full_data = calculate_real_and_gen_proba(pgen_path=clone_pgen_path,
                                 cm_path=cm_path,
                                 run_to_number_of_clonotypes=pd.read_csv(run_to_number_of_clones_path),
                                )
    good_clones = full_data[(full_data['preal'] > preal_threshold) & (full_data['pgen'] > pgen_threshold)].clone

    prepare_clonotype_matrix(cm_path)[['run'] + list(good_clones)].to_csv(save_path)
