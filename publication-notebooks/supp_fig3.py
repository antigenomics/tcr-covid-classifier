#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
# os.chdir('../')
sys.path.append(os.getcwd())

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA
from scipy.stats import fisher_exact, chi2_contingency
from sklearn.manifold import TSNE, MDS
from tqdm import tqdm
from xgboost import XGBClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import f1_score, precision_score, recall_score, roc_curve
from multipy.fwer import sidak, hochberg
import random

from utils.viz_utils import plot_usage_matrix_pca, plot_boxplots_for_usage_matrix, _plot_feature_importance, plot_v_usage_hist, \
                                plot_cluster_map, significant_clones_distribution, plot_results_for_hla_class, plot_generated_to_real_distribution, \
                            plot_olga_cleanup_data, plot_clusters_of_clonotypes, plot_cooccurence_heatmap_with_epitopes_labeling, plot_waterfall_by_column, \
                            plot_feature_importances, plot_volcano
from utils.ml_utils import get_parameters, prepare_data, evaluate_models, split_data_by_batch, split_data, cross_validation_between_batches, make_hla_predictor
from utils.data_utils import prepare_run_column
from utils.stats_utils import evaluate_anova_testing, evaluate_mannwhitneyu_testing, get_top_changed_clonotypes
from utils.clustering_utils import seqs2hamming
from utils.weblogo_utils import create_web_logo
from source.alpha_beta_paired_clones_search import make_metaclone_cm

import warnings
warnings.filterwarnings('ignore')


import importlib
imported_module = importlib.import_module("utils.viz_utils")
importlib.reload(imported_module)


# # COVID-19 associated VDJDb clones selection 

# In[2]:


vdjdb = pd.read_csv('data/vdjdb.txt', sep='\t')


# In[3]:


covid_vdjdb = vdjdb[vdjdb['antigen.species'] == 'SARS-CoV-2']


# In[4]:


covid_vdjdb.head()


# In[5]:


covid_vdjdb[covid_vdjdb.gene == 'TRB']['mhc.a'].apply(lambda x: x.split('-')[-1].split(':')[0]).value_counts()


# In[6]:


covid_vdjdb['mhc.b'].apply(lambda x: x.split('-')[-1].split(':')[0]).value_counts()


# In[7]:


vdjdb_runs = set()
for allele in covid_vdjdb['mhc.a'].apply(lambda x: x.split('-')[-1].split(':')[0]).unique():
    if len(covid_vdjdb[covid_vdjdb['mhc.a'].str.contains(allele)]) > 100:
        try:
            vdjdb_runs = vdjdb_runs.union(set(prepare_run_column(pd.read_csv(f'data/hla_desc/fmba_desc_hla_{allele}.csv')).run))
        except Exception as e:
            print(f'error processing {allele}')
    else:
        print(f'{allele} is too rare')


# In[8]:


for allele in covid_vdjdb['mhc.b'].apply(lambda x: x.split('-')[-1].split(':')[0]).unique():
    if len(covid_vdjdb[covid_vdjdb['mhc.b'].str.contains(allele)]) > 100:
        try:
            vdjdb_runs = vdjdb_runs.union(set(prepare_run_column(pd.read_csv(f'data/hla_desc/fmba_desc_hla_{allele}.csv')).run))
        except Exception as e:
            print(f'error processing {allele}')
    else:
        print(f'{allele} is too rare')


# In[9]:


len(vdjdb_runs)


# In[10]:


covid_vdjdb['vdjdb.score'].value_counts()


# In[11]:


desc = pd.read_csv('data/desc_fmba_not_nan_hla.csv').drop(columns=['Unnamed: 0'])
desc['covid'] = desc.COVID_status.apply(lambda x: 'covid' if x == 'COVID' else 'healthy')
desc


# In[12]:


train_runs = prepare_run_column(desc[~desc.folder.str.lower().str.contains('novaseq6')]).run
test_runs = prepare_run_column(desc[desc.folder.str.lower().str.contains('novaseq6')]).run

train_runs = train_runs[train_runs.isin(vdjdb_runs)]
test_runs = test_runs[test_runs.isin(vdjdb_runs)]


# In[13]:


train_runs.shape


# In[14]:


test_runs.shape


# # VDJDb clones analysis 

# In[15]:


beta_cm = pd.read_csv('data/significant_clone_matrix_fisher_fmba_TRB_vdjdb.csv').drop(columns=['Unnamed: 0'])


# In[16]:


beta_cm


# In[17]:


beta_cm.columns = pd.Series(beta_cm.columns).apply(lambda x: x.split('.')[0])


# In[18]:


beta_cm = beta_cm.loc[:,~beta_cm.columns.duplicated()].copy()


# In[19]:


beta_vdjdb_clones = pd.read_csv('data/significant_clone_matrix_fisher_fmba_TRB_vdjdb.csv').drop(columns=['Unnamed: 0', 'run']).columns


# In[20]:


beta_cm = prepare_run_column(beta_cm)


# In[21]:


significant_clones_distribution(significant_clonotype_matrix=beta_cm[beta_cm.run.isin(vdjdb_runs)], 
                                run_to_number_of_clones=pd.read_csv('data/run_to_number_of_clones_fmba_TRB.csv'), 
                                desc=desc, 
                                data_description='TCRβ biomarkers',
                                by='covid',)


# In[22]:


beta_cm[beta_cm.run.isin(train_runs)]


# In[23]:


fold_change_vdjdb_beta = get_top_changed_clonotypes(clonotype_matrix=beta_cm[beta_cm.run.isin(train_runs)],
                           desc=desc, 
                           pvals=pd.read_csv('data/covid_significant_clone_pvals_fmba_TRB_vdjdb.csv').drop(columns=['Unnamed: 0']), 
                           run_to_number_of_clones=pd.read_csv('data/run_to_number_of_clones_fmba_TRB.csv'),
                           healthy_col='covid', 
                           healthy_label='healthy')
fold_change_vdjdb_beta.dropna()


# In[278]:


plot_volcano(fold_change_vdjdb_beta, pval_threshold=0.2, fold_change_threshold=1.25)


# In[29]:


beta_cm_selected = beta_cm[['run'] + list(
    fold_change_vdjdb_beta[(fold_change_vdjdb_beta.log_fold_change > 1.25) & (fold_change_vdjdb_beta.pval < 0.2)].clone)]


# In[30]:


significant_clones_distribution(significant_clonotype_matrix=beta_cm_selected[beta_cm_selected.run.isin(test_runs)], 
                                run_to_number_of_clones=pd.read_csv('data/run_to_number_of_clones_fmba_TRB.csv'), 
                                desc=desc, 
                                data_description='TCRβ biomarkers',
                                by='covid',)


# In[31]:


beta_covid_vdjdb = vdjdb[(vdjdb['antigen.species'] == 'SARS-CoV-2') & (vdjdb.gene == 'TRB')]\
                        [['cdr3', 'antigen.epitope']].rename(columns={'antigen.epitope': 'cluster'})
beta_covid_vdjdb


# In[32]:


pvals_data = pd.read_csv('data/covid_significant_clone_pvals_fmba_TRB_vdjdb.csv').drop(columns=['Unnamed: 0']).rename(columns={'clone': 'cdr3'})


# In[33]:


pvals_data


# In[34]:


cluster_pvals = beta_covid_vdjdb.merge(pvals_data).groupby(by='cluster', as_index=False).mean().rename(columns={'cluster': 'clone'})
cluster_pvals.clone = 'cluster_' + cluster_pvals.clone
cluster_pvals


# In[35]:


metaclone_beta_cm = make_metaclone_cm(beta_cm, beta_covid_vdjdb[beta_covid_vdjdb.cdr3.isin(beta_cm.columns)])


# In[36]:


metaclone_beta_cm.to_csv('data/clone_matrix_covid_fmba_TRB_metaclone_vdjdb.csv')


# In[37]:


metaclone_beta_cm


# In[38]:


fold_change_vdjdb_beta = get_top_changed_clonotypes(clonotype_matrix=metaclone_beta_cm[metaclone_beta_cm.run.isin(train_runs)],
                           desc=desc, 
                           pvals=cluster_pvals, 
                           log_fold_change_threshold=2, 
                           logp_threshold=1,
                           healthy_col='covid', 
                           healthy_label='healthy')
fold_change_vdjdb_beta.dropna()


# In[39]:


data_with_epi_allele = fold_change_vdjdb_beta.dropna().sort_values(by='log_fold_change')
data_with_epi_allele.clone = data_with_epi_allele.clone.apply(lambda x: x.split('_')[1])
data_with_epi_allele = data_with_epi_allele.rename(columns={'clone': 'antigen.epitope'}).merge(vdjdb[['antigen.epitope', 'mhc.a']])
data_with_epi_allele[data_with_epi_allele.log_fold_change > 1].drop_duplicates()['mhc.a'].value_counts()


# In[40]:


data_with_epi_allele[data_with_epi_allele.log_fold_change < 100]


# In[41]:


plot_volcano(fold_change_vdjdb_beta, pval_threshold=1, fold_change_threshold=1)


# In[42]:


beta_cm_selected_meta = metaclone_beta_cm[['run'] + list(
    fold_change_vdjdb_beta[(fold_change_vdjdb_beta.log_fold_change > 1) & (fold_change_vdjdb_beta.pval <= 1)].clone)]


# In[43]:


beta_cm_selected_meta


# In[44]:


metaclone_beta_cm['cluster_AFLLFLVLI']


# In[45]:


alpha_cm = pd.read_csv('data/significant_clone_matrix_fisher_fmba_TRA_vdjdb.csv').drop(columns=['Unnamed: 0'])
alpha_vdjdb_clones = pd.read_csv('data/significant_clone_matrix_fisher_fmba_TRA_vdjdb.csv').drop(columns=['Unnamed: 0', 'run']).columns
alpha_covid_vdjdb = vdjdb[(vdjdb['antigen.species'] == 'SARS-CoV-2') & (vdjdb.gene == 'TRA')]\
                        [['cdr3', 'antigen.epitope']].rename(columns={'antigen.epitope': 'cluster'})


# In[46]:


alpha_cm = prepare_run_column(alpha_cm)


# In[47]:


significant_clones_distribution(significant_clonotype_matrix=alpha_cm[alpha_cm.run.isin(vdjdb_runs)], 
                                run_to_number_of_clones=pd.read_csv('data/run_to_number_of_clones_fmba_TRA.csv'), 
                                desc=desc, 
                                data_description='TCRα biomarkers',
                                by='covid',)


# In[48]:


fold_change_vdjdb_alpha = get_top_changed_clonotypes(clonotype_matrix=alpha_cm[alpha_cm.run.isin(train_runs)],
                           desc=desc, 
                           pvals=pd.read_csv('data/covid_significant_clone_pvals_fmba_TRA_vdjdb.csv').drop(columns=['Unnamed: 0']), 
                           log_fold_change_threshold=2, 
                           logp_threshold=1,
                           healthy_col='covid', 
                           healthy_label='healthy')
fold_change_vdjdb_alpha.dropna()


# In[49]:


plot_volcano(fold_change_vdjdb_alpha, pval_threshold=0.2, fold_change_threshold=1.4)


# In[50]:


alpha_cm_selected = alpha_cm[['run'] + list(
    fold_change_vdjdb_alpha[(fold_change_vdjdb_alpha.log_fold_change > 1.4) & (fold_change_vdjdb_alpha.pval < 0.2)].clone)]


# In[51]:


significant_clones_distribution(significant_clonotype_matrix=alpha_cm_selected[alpha_cm_selected.run.isin(test_runs)], 
                                run_to_number_of_clones=pd.read_csv('data/run_to_number_of_clones_fmba_TRA.csv'), 
                                desc=desc, 
                                data_description='TCRα biomarkers',
                                by='covid',)


# In[52]:


metaclone_alpha_cm = make_metaclone_cm(alpha_cm_selected, alpha_covid_vdjdb[alpha_covid_vdjdb.cdr3.isin(alpha_cm_selected.columns)])
metaclone_alpha_cm.to_csv('data/clone_matrix_covid_fmba_TRA_metaclone_vdjdb.csv')


# # Training separate α/β classifiers

# In[201]:


runs = prepare_run_column(pd.read_csv('data/standardized_usage_matrix_fmba_TRA.csv'))


# In[205]:


runs = runs[runs.run.isin(vdjdb_runs)]


# In[208]:


data_beta = prepare_data(run_to_number_of_clones_path='data/run_to_number_of_clones_fmba_TRB.csv',
                     desc_path='data/standardized_usage_matrix_fmba_TRA.csv',
                     clonotype_matrix_path='data/significant_clone_matrix_fisher_fmba_TRB_vdjdb.csv',
                     hla_keys_path='data/hla_keys.csv',
                     make_freq=True, 
                     make_all_features_bool=False,
                     use_hla_clones=False, 
                     use_hla_bool=True, 
                     use_standardization=True,
                     raw_target_column='covid',
                     raw_target_clumn_success_label='covid',).loc[runs.index]


# In[209]:


data_beta


# In[210]:


data_beta = data_beta[list(beta_cm_selected.columns[1:]) + ['covid', 'is_test_run']]


# In[211]:


data_beta.shape


# In[212]:


X_train, y_train, X_test, y_test = split_data(data=data_beta, y_column='covid', by='is_test_run')
best_clfs_beta = evaluate_models(X_train, y_train, X_test, y_test, get_parameters(), scoring_function='f1_weighted', debug=True)


# In[213]:


best_clfs_beta[0]['svm'].predict(X_test)


# In[214]:


pred = pd.DataFrame({'test': y_test, 'pred':best_clfs_beta[0]['svm'].predict(X_test)})
pred['correct'] = pred.test == pred.pred


# In[215]:


pred.correct.value_counts()


# In[216]:


best_clfs_beta[0]['svm'].score(X_test, y_test)


# In[217]:


proba_labels = y_test.apply(lambda x: 'healthy' if x == 0 else 'covid')
probability_df = pd.DataFrame({'beta_proba': best_clfs_beta[0]['svm'].predict_proba(X_test)[::, 1], 'covid': proba_labels})
plot_waterfall_by_column(probability_df, proba_column='beta_proba', label_column='covid')


# In[219]:


data_alpha = prepare_data(run_to_number_of_clones_path='data/run_to_number_of_clones_fmba_TRA.csv',
                     desc_path='data/standardized_usage_matrix_fmba_TRA.csv',
                     clonotype_matrix_path='data/significant_clone_matrix_fisher_fmba_TRA_vdjdb.csv',
                     hla_keys_path='data/hla_keys.csv',
                     make_freq=True,
                     use_hla_clones=False,
                     use_hla_bool=True,
                     use_standardization=True,
                     raw_target_column='covid',
                     raw_target_clumn_success_label='covid',).loc[runs.index]


# In[220]:


data_alpha = data_alpha[list(alpha_cm_selected.columns[1:]) + ['covid', 'is_test_run']]


# In[221]:


X_train_alpha, y_train_alpha, X_test_alpha, y_test_alpha = split_data(data=data_alpha, y_column='covid', by='is_test_run')

best_clfs_alpha = evaluate_models(X_train_alpha, y_train_alpha, 
                                  X_test_alpha, y_test_alpha, 
                                  get_parameters(), scoring_function='f1_weighted', debug=True)


# # Training joint αβ classifier

# In[222]:


data_alpha_for_joint = prepare_data(run_to_number_of_clones_path='data/run_to_number_of_clones_fmba_TRA.csv',
                     desc_path='data/standardized_usage_matrix_fmba_TRA.csv',
                     clonotype_matrix_path='data/significant_clone_matrix_fisher_fmba_TRA_vdjdb.csv',
                     hla_keys_path='data/hla_keys.csv',
                     make_freq=True,
                     use_hla_clones=False,
                     use_hla_bool=False,
                     use_standardization=True,
                     raw_target_column='covid',
                     raw_target_clumn_success_label='covid',
                     metadata_columns=['project', 'is_test_run']).loc[runs.index]
print('alpha ready')
data_beta_for_joint = prepare_data(run_to_number_of_clones_path='data/run_to_number_of_clones_fmba_TRB.csv',
                     desc_path='data/standardized_usage_matrix_fmba_TRA.csv',
                     clonotype_matrix_path='data/significant_clone_matrix_fisher_fmba_TRB_vdjdb.csv',
                     hla_keys_path='data/hla_keys.csv',
                     make_freq=True, 
                     use_hla_clones=False, 
                     use_hla_bool=True, 
                     use_standardization=True,
                     raw_target_column='covid',
                     raw_target_clumn_success_label='covid',).loc[runs.index]


# In[223]:


data_joint = pd.concat([data_beta_for_joint, data_alpha_for_joint.drop(columns=['covid', 'is_test_run'])], axis=1)


# In[224]:


data_joint = data_joint[list(alpha_cm_selected.columns[1:]) + list(beta_cm_selected.columns[1:]) + ['covid', 'project', 'is_test_run']]


# In[225]:


data_joint


# In[226]:


X_train_joint, y_train_joint, X_test_joint, y_test_joint = split_data(data=data_joint.drop(columns=['project']), y_column='covid', by='is_test_run')

best_clfs_joint = evaluate_models(X_train_joint, y_train_joint, X_test_joint, y_test_joint, 
                                  get_parameters(), scoring_function='f1_weighted', debug=True)


# In[227]:


all_clfs = {
     'beta': best_clfs_beta,
     'alpha': best_clfs_alpha,
     'joint': best_clfs_joint
}
X_test_data={
     'beta': X_test,
     'alpha': X_test_alpha,
     'joint': X_test_joint
}
y_test_data={
     'beta': y_test,
     'alpha': y_test_alpha,
     'joint': y_test_joint
}
X_train_data={
     'beta': X_train,
     'alpha': X_train_alpha,
     'joint': X_train_joint
}
y_train_data={
     'beta': y_train,
     'alpha': y_train_alpha,
     'joint': y_train_joint
}


# In[228]:


data_joint


# In[229]:


model_df = []
f1_df = []
data_type_df = []
for data_type in ['beta','alpha','joint']:
    for model in ['svm', 'ab', 'knn', 'rfc', 'mlpclassifier', 'xgboost']:
        cur_score = f1_score(y_test_data[data_type], all_clfs[data_type][0][model].predict(X_test_data[data_type]))
        model_df.append(model)
        data_type_df.append(data_type)
        f1_df.append(cur_score)
comparison_df = pd.DataFrame({'model':model_df, 'f1': f1_df, 'data_type':data_type_df})


# In[230]:


comparison_df


# In[231]:


fig, ax = plt.subplots()
sns.catplot(data=comparison_df, x="model", y="f1", hue="data_type", kind="point", ax=ax)
plt.close(1)


# In[232]:


sns.stripplot(data=comparison_df, x="model", y="f1", hue="data_type",)


# In[233]:


all_clfs['joint']


# # Making metrics dataframe

# In[253]:


clf_name = []
clf_type = []
f1 = []
precision = []
recall = []
for key, clfs in all_clfs.items():
    clf_name.append(key)
    clf_type.append(clfs[2])
    best_clf = clfs[0]['xgboost']
    clf_predictions = best_clf.predict(X_test_data[key])
    f1.append(f1_score(y_test_data[key], clf_predictions))
    precision.append(precision_score(y_test_data[key], clf_predictions))
    recall.append(recall_score(y_test_data[key], clf_predictions))


# In[254]:


plotting_df = pd.DataFrame({
    'classifier': clf_name, 
    'best_classifier_predictor': clf_type,
    'f1_score': f1,
    'precision_score': precision,
    'recall_score': recall
})


# In[255]:


plotting_df


# In[256]:


plotting_df = plotting_df.applymap(lambda x: round(x, 2) if isinstance(x, float) else x)


# # Preparing data for proba comparison

# In[257]:


beta_predictions = all_clfs['beta'][0]['xgboost'].predict_proba(X_test_data['beta'])[::,1]
alpha_predictions = all_clfs['alpha'][0]['xgboost'].predict_proba(X_test_data['alpha'])[::,1]

proba_labels = y_test.apply(lambda x: 'healthy' if x == 0 else 'covid')
probability_df = pd.DataFrame({
    'beta_proba': beta_predictions,
    'alpha_proba': alpha_predictions,
    'covid': proba_labels
})


# In[258]:


probability_df


# # One folder out CV

# In[259]:


data_joint


# In[260]:


desc_for_projects = pd.read_csv('data/desc_fmba_new_split.csv').drop(columns=['Unnamed: 0'])[['run', 'folder']]
desc_for_projects.head()


# In[261]:


data_joint_proj = data_joint.copy()


# In[262]:


data_joint_proj['project'] = desc_for_projects.folder.apply(lambda x: x.replace('_DNA', '').split('/')[-1].split('_')[-1])


# In[263]:


data_joint_proj.project.unique()


# In[264]:


data_joint_proj.head()


# In[265]:


data_joint_proj


# In[266]:


metrics_df = []
for metrics, metrics_name in zip([f1_score, precision_score, recall_score], ['f1', 'precision', 'recall']):
    metrics_df.append(cross_validation_between_batches(best_clfs_joint[0]['xgboost'], 
                                 data_joint_proj.drop(columns=['is_test_run']), 
                                 [x for x in data_joint_proj.project.unique()], 
                                 y_column='covid', 
                                 batch_column='project', 
                                 metrics=metrics, 
                                 metrics_name=metrics_name, 
                                 return_metrics_results=True,
                                 debug=True
                                 ))


# In[267]:


def make_score_column(df):
    df['score'] = df[df.columns[1]]
    df['metrics'] = df.columns[1]
    df = df.drop(columns=[df.columns[1]])
    return df

metrics_df = pd.concat([make_score_column(metrics_df[i]) for i in range(3)])
metrics_df.folder = metrics_df.folder.apply(lambda x: x.replace('_DNA', '').split('/')[-1].split('_')[-1])


# # Plotting

# In[268]:


for key in all_clfs:
    all_clfs[key][0]['xgboost'].fit(X_train_data[key], y_train_data[key])


# In[269]:


all_clfs.keys()


# In[279]:


fig = plt.figure(figsize=(15, 20))
gs = GridSpec(nrows=9, 
              ncols=3)
font_size=20
delta_x=-0.1
delta_y=1.14

########################################################

ax0 = fig.add_subplot(gs[:2, 0])
plot_volcano(fold_change_vdjdb_beta, pval_threshold=0.2, fold_change_threshold=1.25, ax=ax0)

ax1 = fig.add_subplot(gs[:2, 1])
ax2 = fig.add_subplot(gs[:2, 2])

significant_clones_distribution(significant_clonotype_matrix=beta_cm_selected[beta_cm_selected.run.isin(test_runs)], 
                                run_to_number_of_clones=pd.read_csv('data/run_to_number_of_clones_fmba_TRB.csv'), 
                                desc=desc, 
                                data_description='TCRβ biomarkers',
                                by='covid',
                                reg_ax=ax1, hist_ax=ax2)

ax0.text(delta_x, delta_y, 'A',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax0.transAxes,
        size=font_size)
ax1.text(delta_x, delta_y, 'B',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax1.transAxes,
        size=font_size)
ax2.text(delta_x, delta_y, 'C',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax2.transAxes,
        size=font_size)

########################################################

ax3 = fig.add_subplot(gs[2:4, 0])
plot_volcano(fold_change_vdjdb_alpha, pval_threshold=0.2, fold_change_threshold=1.4, ax=ax3)

ax4 = fig.add_subplot(gs[2:4, 1])
ax5 = fig.add_subplot(gs[2:4, 2])

significant_clones_distribution(significant_clonotype_matrix=alpha_cm_selected[alpha_cm_selected.run.isin(test_runs)], 
                                run_to_number_of_clones=pd.read_csv('data/run_to_number_of_clones_fmba_TRA.csv'), 
                                desc=desc, 
                                data_description='TCRα biomarkers',
                                by='covid',
                                reg_ax=ax4, hist_ax=ax5)

ax3.text(delta_x, delta_y, 'D',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax3.transAxes,
        size=font_size)
ax4.text(delta_x, delta_y, 'E',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax4.transAxes,
        size=font_size)
ax5.text(delta_x, delta_y, 'F',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax5.transAxes,
        size=font_size)

########################################################

for i, (letter, metrics) in enumerate(zip(['J', None, None], ['f1_score', 'precision_score', 'recall_score'])):
    
    ax = fig.add_subplot(gs[6, i])
    plotting_df[['classifier', metrics]].set_index('classifier').plot(kind='bar', ax=ax)
    ax.bar_label(ax.containers[0])
    ax.tick_params(labelrotation=20)
    ax.get_legend().remove()
    ax.set_ylabel(metrics)
    ax.set_title(f"{metrics.replace('_', ' ')}")
    ax.tick_params(labelrotation=20)
    ax.set_ylim(0, 1)
    if letter is not None:
        ax.text(delta_x, delta_y, letter,
             horizontalalignment='left',
             verticalalignment='top',
            transform=ax.transAxes,
                size=font_size)

########################################################
ax6 = fig.add_subplot(gs[4:6, 1])

for key in all_clfs:
    PrecisionRecallDisplay.from_estimator(
        all_clfs[key][0][all_clfs[key][2]], X_test_data[key], y_test_data[key], name=key, ax=ax6
    )
    # fpr, tpr, _ = roc_curve(y_test_data[key], all_clfs[key][0][all_clfs[key][2]].predict_proba(X_test_data[key])[::,1])
    # ax3.plot(fpr,tpr, label=key)

ax6.set_ylabel('Precision')
ax6.set_xlabel('Recall')
ax6.plot([0, 1], [1, 1], linestyle='dashed', color='grey')
ax6.plot([1, 1], [1, 0], linestyle='dashed', color='grey')
ax6.legend()
# ax3.set_title('ROC-curve')
ax6.text(delta_x, delta_y, 'H',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax6.transAxes,
        size=font_size)

########################################################
ax7 = fig.add_subplot(gs[4:6, 2])

sns.stripplot(data=metrics_df, x="metrics", y="score", hue="folder", s=10, ax=ax7)
sns.boxplot(data=metrics_df, x="metrics", y="score", color='azure', ax=ax7)
# ax7.set_title('Scores across batches for α+β metaclone classifier')
ax7.set_ylim(0.5, 1.1)
ax7.text(delta_x, delta_y, 'I',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax7.transAxes,
        size=font_size)

########################################################

ax9 = fig.add_subplot(gs[4:6, 0])
sns.stripplot(data=comparison_df, x="model", y="f1", hue="data_type", ax=ax9, s=10)
sns.boxplot(data=comparison_df, x="model", y="f1", ax=ax9, color='azure')
ax9.text(delta_x+0.08, delta_y+0.05, 'G',
     horizontalalignment='left',
     verticalalignment='top',
    transform=ax9.transAxes,
        size=font_size)
# ax9.set_title('Comparison of f1-score across models')
ax9.tick_params(labelrotation=20)
ax9.set_ylim(0.3, 0.8)
ax9.legend(ncol=3,  loc='upper center')

########################################################

plt.tight_layout()
plt.savefig("figures/supp_fig3.svg")
plt.savefig("figures/supp_fig3.pdf")
plt.savefig("figures/supp_fig3.png")
plt.show()


# Figure 4. Analysis of machine learning approaches applied to FMBA TCRβ and TCRα biomarkers. 
# 
# A, B, C. Distribution of target metrics (f1-score, precision, recall) for all the evaluated models.
# 
# D. ROC-curve plot for all the evaluated models.
# 
# E. The sctterplot describing the relation of the probabilities to label the sample as the COVID-19 positive for TRCα and TCRβ based classifiers.
# 
# F. The sctterplot describing the relation of the probabilities to label the sample as the COVID-19 positive for classifier based on both TRCα,TCRβ biomarkers and TRCα,TCRβ metaclone cluster features.
# 
# G. The waterfall plot representing the probability of each sample to be labeld as COVID-19 positiove (> 0) or healthy (< 0). Samples coming from healthy donors are colored with blue, COVID-19 samples are colored with orange.
# 
# H. Evaluation of target metrics (f1-score, precision, recall) for one batch out cross validation.
# 
# I. Feature importance plot for the XGBoost classifier based on TRCα and TCRβ based biomarkers and HLA presence features.

# In[ ]:





# In[ ]:




