import pandas as pd
from sklearn.model_selection import train_test_split

if __name__ == "__main__":
    um_path = snakemake.input[0]
    um = pd.read_csv(um_path).drop(columns=['Unnamed: 0'])
    if snakemake.params.test_selection_method == 'batch':
        um['is_test_run'] = um.project.apply(lambda x: snakemake.params.test_batch.lower() in x.lower())
        um.to_csv(snakemake.output[0])
    elif snakemake.params.test_selection_method == 'percent':
        healthy_um = um[um.covid == 'healthy'].sort_values(by='run')
        covid_um = um[um.covid != 'healthy'].sort_values(by='run')
        healthy_train_um, healthy_test_um = train_test_split(healthy_um, test_size=snakemake.params.test_percent, shuffle=False)
        covid_train_um, covid_test_um = train_test_split(covid_um, test_size=snakemake.params.test_percent, shuffle=False)
        train_um = pd.concat([healthy_train_um, covid_train_um])
        test_um = pd.concat([healthy_test_um, covid_test_um])
        train_um['is_test_run'] = False
        test_um['is_test_run'] = True
        pd.concat([train_um, test_um]).to_csv(snakemake.output[0])
