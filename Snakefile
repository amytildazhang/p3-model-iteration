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
feature_types = [Path(x).parent.name for x in glob.glob('data/raw/features/*/*.tsv')]

# drug response input files
response_files = os.listdir('data/raw/response')

#
# Rules
#
rule create_training_sets:
    input: expand("data/raw/response/{response}", response=response_files),
           expand("data/processed/features/{feature_type}/{feature}_pca.tsv", zip, feature_type=feature_types, feature=features)

#
# Load data
#
rule pca_project:
    input: 'data/raw/features/{feature_type}/{feature}.tsv'
    output: touch('data/processed/features/{feature_type}/{feature}_pca.tsv')
    #script: 'scripts/pca_project.R'
