#
# Multiple Myeloma Dose Response Prediction Model Iteration
# V. Keith Hughitt
# June 2019
#
import glob
from pathlib import Path

configfile: "config.yaml"

workdir: "/data/projects/nih/p3-models"

# expected outputs
features = [Path(x).stem for x in glob.glob('data/raw/features/*/*.tsv')]

# create a list with each feature type repeated based on the number of input datasets
# associated with that feature type (e.g. ['cnv', 'cnv', 'rnaseq', 'variants', ...]
#feature_types = ['cnv' for i in range(len(os.listdir('input/features/cnv')))] +
#                ['rnaseq' for i in range(len(os.listdir('input/features/rnaseq')))] 
#                ['variants' for i in range(len(os.listdir('input/features/variants')))] 
feature_types = [Path(x).parent.name for x in glob.glob('data/raw/features/*/*.tsv')]

drug_datasets = [Path(x).stem for x in glob.glob('data/raw/response/*.tsv')]

#
# Rules
#
rule all:
    input: expand("data/processed/response/{drug_dataset}/{drug}.tsv", drug_dataset=drug_datasets, drug=config['target_drugs']),
           expand("data/processed/features/{feature_type}/{feature}_pca.tsv", zip, feature_type=feature_types, feature=features)

#
# Load data
#
rule split_drugs:
    input:  expand("data/raw/response/{drug_dataset}.tsv", drug_dataset=drug_datasets)
    output: expand("data/processed/response/{drug_dataset}/{drug}.tsv", drug_dataset=drug_datasets, drug=config['target_drugs'])
    script: 'scripts/split_drugs.R'

rule pca_project:
    input: 'data/raw/features/{feature_type}/{feature}.tsv'
    output: 'data/processed/features/{feature_type}/{feature}_pca.tsv'
    script: 'scripts/pca_project.R'
