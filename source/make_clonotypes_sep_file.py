import pandas as pd
import os

os.chdir('/home/ekvlasova/covid-t-cell-classifier')


def make_sep_file(path, save_path):
    df = pd.read_csv(path)[['clone']]
    df.to_csv(save_path, header=None)


if __name__ == "__main__":
    # make_sep_file(path='data/clonotype_matrix_fmba_top_500k_1_mismatch.csv',
    #               save_path='data/olga_pgen/fmba_top_500k_1_mismatch.csv')
    make_sep_file(path='data/alpha/covid_clones_all_hla_top_500k_1_mismatch_fisher.csv',
                  save_path='data/olga_pgen/alpha_covid_clones.csv')
