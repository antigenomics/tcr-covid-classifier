# My TODOs
3. Сделать картинку 3
4. Сделать всю подготовку для адаптив
5. Сделать картинку 5


# TCR alpha and beta based classifier for COVID-19 status prediction
The repository contains code for paper "Predicting SARS-CoV-2 exposure using T-cell 
repertoire sequencing and machine learning".
The overall pipeline is shown below: 
![Pipeline overview](figures/fig0.png)
## Data availability
All the data used in this research was uploaded to Zenodo platform.
**TODO**
## MetaClone classifier
### Prerequisites
We expect that you have an access to a cluster with at least 
* 24 threads 
* 16gb RAM 

You need to have Python, conda and snakemake installed to run the code.
We use OLGA software for TCR generation probability calculations.

### Installation
1. `git clone https://github.com/antigenomics/tcr-covid-classifier.git`
2. `conda env create -f environment.yml`
3. `conda activate tcr-classifier`

### Run the whole pipeline
The code is organized into a snakemake pipeline.
You need to run the following command from the root directory:
```cmd
snakemake --cores 48 
```
You can specify less cores, but the code would be requiring more time to execute.
Using 16gb RAM and 48 CPUs the code is expected to run for ~2 hours.  
