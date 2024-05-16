# TCR alpha and beta based classifier for COVID-19 status prediction
The repository contains code for article "Robust SARS-CoV-2 exposure 
biomarkers in population detected using T-cell repertoire profiling".

The overall pipeline is shown below: 
![Pipeline overview](figures/fig0.png)
The visualization of the biomarkers inferred can be found [here](https://covidbiomarkers.cdr3.net/) (under construction). 
For a more detailed description of the dataset see Zenodo record.
## Data availability
All the data generated in this research was uploaded to [Zenodo platform](https://zenodo.org/records/8362803) (a cohort of >1200 healthy donors and COVID-19 patients).
The dataset contains ~40 millions TCR alpha clonotypes and ~20 millions TCR beta clonotypes altogether. 
This makes our dataset one of the largest repertoire sequencing datasets published to date.
For more information about the dataset see Zenodo record.
## Biomarkers information
The list of the inferred biomarkers, i.e. clonotype with significant association with COVID-19 status, 
can be found [here](figures/covid_clonotypes.csv). Each biomarker has a mapping to the cluster 
(groups of biomarkers with similar CDR3 sequences) it belongs to, clusters are 0-based numbered 
for alpha and beta chains separately.
Each biomarker record contains the following information:
* cdr3 amino acid sequence
* cluster id (0..57 for beta, 0..114 for alpha)
* chain (alpha/beta)
* most frequent v/j segments which occur in the data for this CDR3 sequence (might be several if they occur with equal frequency)

The clustering summary can be found [here](figures/clustering_summary_TRA.csv) for TCR alpha chain and 
[here](figures/clustering_summary_TRB.csv) for TCR beta chain.
The clustering summary files contain the following columns:
* cluster, cluster_size -- info on cluster id and its size
* information on the associated antigen epitope according to VDJdb (empty if associations are not found) 
  * antigen_epitope
  * antigen_species 
* num_epi_associated_clones -- number of clonotypes which are associated with the epitope
* pval -- p-value of the association
* all_clust_count -- number of epitope associations within all clusters
* enrichment_score -- fraction of epitope association occurences in the cluster compared to all the clusters 
* mhc of the epitope in VDJdb
  * mhc_a
  * mhc_b
* fraction_epi_associated_clones -- num_epi_associated_clones divided by cluster size 
* information on cluster incidence in training set, i.e. number of samples from training set containing at least one biomarker from the cluster 
  * num_samples_with_cluster_train 
  * num_samples_with_cluster_train_healthy
  * num_samples_with_cluster_train_covid
* cluster incidence in test set
  * num_samples_with_cluster_test
  * num_samples_with_cluster_test_healthy
  * num_samples_with_cluster_test_covid
* cluster frequency in train set 
  * fraction_samples_with_cluster_train
  * fraction_of_healthy_samples_in_cluster_train
  * fraction_of_covid_samples_in_cluster_train
* cluster frequency in test set
  * fraction_samples_with_cluster_test
  * fraction_of_healthy_samples_in_cluster_test
  * fraction_of_covid_samples_in_cluster_test
## MetaClone classifier
### Prerequisites
We expect that you have an access to a cluster with at least 
* 24 threads 
* 16gb RAM 

You need to have Python, conda and snakemake installed to run the code.
We use [OLGA software](https://github.com/statbiophys/OLGA) for TCR generation probability calculations.

### Installation
1. `git clone https://github.com/antigenomics/tcr-covid-classifier.git`
2. `conda env create -f environment.yml`
3. `conda activate tcr-classifier`
4. `pip install logomaker openpyxl`

### Run the whole pipeline
The code is organized into a snakemake pipeline.
You need to run the following command from the root directory:
```cmd
snakemake --cores 48 
```
You can specify less cores, but the code would be requiring more time to execute.
Using 16gb RAM and 48 CPUs the code is expected to run for ~12 hours.  

### TCR biomarkers search
The biomarkers search is implemented in ``source/tests_analysis/covid_test_250k.py``. 
You will need to prepare the data firstly. Data preparation includes the following steps:
1. Segment usage matrix creation (``source/v_j_usage_matrix_creation.py``)
2. Usage matrix standardization if data was prepared in more than one batch (`source/usage_matrix_preprocessing.py`)
3. Data resampling (`source/downsampling_procedure`) which is required in case when step 2 was used
4. Public clones extraction and clone matrix creation (`source/clonotypes_extraction_procedure.py` and `source/clonotype_matrix_creation.py`)

### Evaluation of found biomarkers
The biomarkers are analyzed on how they can be clustered together within one amino acid change. In addition to that we 
implemented a pipeline to search for associations of the found clusters in VDJdb. 
This analysis is performed in ``publication-notebooks/fig2.ipynb``.

You can try the analysis yourself by browsing our [website](https://covidbiomarkers.cdr3.net/).

### Classification
The classification pipeline includes cross-validation between batches of data and comparing different models.
The more precise description of the analysis performed can be found in ``publication-notebooks/fig4.ipynb`` and 
`publication-notebooks/fig5.ipynb` (for classification between data derived with different platforms).

### HLA specific biomarkers search increases classification quality
The analysis of A*02 classification quality vs all HLA data can be found in ``publication-notebooks/fig6.ipynb``. 
The analysis of rare HLA pattern found for DRB1\*16 and DQB1\*05 positive donors is located in ``publication-notebooks/fig4.ipynb``
