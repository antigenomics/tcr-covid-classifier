from multiprocessing import Manager, Pool

import pandas as pd


def process_all_files(desc_df, raw_data_path, threads=32, make_read_col=False):
    global process_one_file
    run_to_clone_number = Manager().dict()
    run_to_read_number = Manager().dict()

    def process_one_file(run):
        try:
            file = f'{raw_data_path}/{run}'
            cur_data = pd.read_csv(file, sep=',' if file.split('.')[-1] == 'csv' else '\t')
            run_to_clone_number[run] = cur_data.shape[0]
            run_to_read_number[run] = cur_data['count'].sum()
        except Exception as e:
            pass

    if 'run' in desc_df.columns:
        runs = list(desc_df['run'])
    else:
        runs = list(desc_df['file_name'])
    print(len(runs))
    with Pool(threads) as p:
        p.map(process_one_file, runs)
    data_dict = {'run': [], 'number_of_clones': []}
    if make_read_col:
        data_dict['number_of_reads'] = []
    for x, y in run_to_clone_number.items():
        data_dict['run'].append(x)
        data_dict['number_of_clones'].append(y)
        if make_read_col:
            data_dict['number_of_reads'].append(run_to_read_number[x])
    return pd.DataFrame(data=data_dict)


if __name__ == "__main__":
    file = snakemake.input[0]
    um = pd.read_csv(file, sep=',' if file.split('.')[-1] == 'csv' else '\t')
    print(um)
    process_all_files(desc_df=um,
                      raw_data_path=snakemake.input[1],
                      threads=snakemake.threads,
                      make_read_col=snakemake.params.make_read_count_col).to_csv(snakemake.output[0], index=False)

