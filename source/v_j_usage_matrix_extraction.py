import os
from collections import defaultdict

import pandas as pd
from tqdm import tqdm


def generate_v_gene_usage_vector(sample):
    sample['v_gene'] = sample.v
    sample['j_gene'] = sample.j
    v_gene_usage_df = sample[['v_gene', 'j_gene', 'cdr3nt']].sort_values(by='v_gene').groupby(['v_gene', 'j_gene'],
                                                                                              as_index=False).count()
    vj_gene_to_usage_dict = {}
    for v, j, count in zip(v_gene_usage_df['v_gene'], v_gene_usage_df['j_gene'], v_gene_usage_df['cdr3nt']):
        vj_gene_to_usage_dict[v, j] = count
    return vj_gene_to_usage_dict


def create_joint_vj_usage_matrix(repertoires_list_project, runs_list_project, output):
    vs = []
    js = []
    runs = []
    projects = []
    counts = []
    for project in repertoires_list_project:
        for i, run_repertoire in enumerate(repertoires_list_project[project]):
            run = runs_list_project[project][i]
            for (v, j), count in run_repertoire.items():
                vs.append(v)
                js.append(j)
                runs.append(run)
                projects.append(project)
                counts.append(count)
    pd.DataFrame({'v': vs,
                  'j': js,
                  'count': counts,
                  'run': runs,
                  'project': projects}).to_csv(output)


def process_all_files(all_runs_path,
                      output_path='/home/ekvlasova/covid-t-cell-classifier/data',
                      raw_data_path='/projects/fmba_covid/1_data_links/fmba_2021',
                      gene='TRB',
                      platform='fmba',
                      create_sep_files=False):
    print(os.getcwd())
    runs = pd.read_csv(all_runs_path, sep='\t')
    runs = runs[runs['file.name'].str.contains(gene)]
    repertoires_list_project = defaultdict(list)
    runs_list_project = defaultdict(list)
    for run, project in tqdm(zip(runs['file.name'], runs['folder']), total=len(runs),
                             desc='Extracting v gene usage vectors'):
        data = pd.read_csv(f'{raw_data_path}/{run}')
        run_repertoire = generate_v_gene_usage_vector(data)
        repertoires_list_project[project].append(run_repertoire)
        runs_list_project[project].append(run)

    create_joint_vj_usage_matrix(repertoires_list_project, runs_list_project,
                                 output=output_path + f'/vj_usage_matrix_{platform}_{gene}.csv')

    if create_sep_files:
        for project in repertoires_list_project:
            vj_to_count = defaultdict(int)
            for run_repertoire in repertoires_list_project[project]:
                for (v, j), count in run_repertoire.items():
                    vj_to_count[v, j] += count
            vs = [v for (v, j) in vj_to_count]
            js = [j for (v, j) in vj_to_count]
            counts = [vj_to_count[k] / len(repertoires_list_project[project]) for k in vj_to_count]
            print(project, f'{output_path}_{project}.csv')
            pd.DataFrame({'v': vs, 'j': js, 'count': counts}).to_csv(f'{output_path}_{project.replace("/", "")}.csv')


if __name__ == "__main__":
    if 'snakemake' in globals():
        print('Snakemake launch')
        process_all_files(all_runs_path=snakemake.params['all_runs_path'],
                          output_path=snakemake.params['output_path'],
                          raw_data_path=snakemake.params['raw_data_path'],
                          gene=snakemake.params['gene'],
                          platform=snakemake.params['platform'])
    else:
        print('Default launch')
        process_all_files(all_runs_path='/projects/fmba_covid/1_data_links/fmba_2021.txt')
