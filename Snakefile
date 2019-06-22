#
# Multiple Myeloma Dose Response Prediction Model Iteration
# V. Keith Hughitt
# June 2019
#
import glob
import os
from pathlib import Path

configfile: "config.yaml"

# expected outputs
data_dir = os.path.join('data', config['name'])

feature_files = glob.glob(os.path.join(data_dir, 'raw/features/*/*.tsv.gz'))
features = [Path(x).stem.replace('.tsv', '') for x in feature_files]

# create a list with each feature type repeated based on the number of input datasets
# associated with that feature type (e.g. ['cnv', 'cnv', 'rnaseq', 'variants', ...]
feature_types = [Path(x).parent.name for x in feature_files]

# drug response input files
response_files = os.listdir(os.path.join(data_dir, 'raw/response'))

#
# Rules
#
rule create_training_sets:
    input: expand(os.path.join(data_dir, "raw/response/{response}"), response=response_files),
           expand(os.path.join(data_dir, "processed/features/{feature_type}/{feature}_pca.tsv.gz"), zip, feature_type=feature_types, feature=features)

#
# Load data
#
rule pca_project:
    input: os.path.join(data_dir, 'raw/features/{feature_type}/{feature}.tsv.gz')
    output: os.path.join(data_dir, 'processed/features/{feature_type}/{feature}_pca.tsv.gz')
    script: 'scripts/pca_project.R'
