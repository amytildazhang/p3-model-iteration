#
# Dose Response Prediction Model Iteration
# V. Keith Hughitt, Amy Zhang
# June 2019
#
# Pipeline steps:
#
# 1. Split data into CV folds
# 2. Gene set projection [Optional]
# 3. Early dimension reduction [Optional]
# 4. Feature filtering (unsupervised)
# 5. Training set construction
# 6. Late dimension reduction [Optional]
# 7. Feature selection (supervised)
# 8. Model training 
#
import glob
import os
import pandas as pd
from os.path import join
from pathlib import Path
from sklearn.model_selection import RepeatedKFold

# report debug mode status
if config['debug']:
    print("=== DEBUG MODE ENABLED ===")

# base data directory
input_dir = join('data', config['name'])
output_dir = join('output', config['name'], str(config['version']))

# load sample metadata
samples = pd.read_csv(join(input_dir, 'metadata', 'samples.tsv'), sep='\t')

# drug response input files
response_files = [Path(x).name for x in glob.glob(join(input_dir, 'response/*.tsv.gz'))]
drug_names = [x.replace('.tsv.gz', '') for x in response_files] 

#
# Cross-validation setup
#
rkf = RepeatedKFold(n_splits=config['cross_validation']['num_splits'],
                    n_repeats=config['cross_validation']['num_repeats'],
                    random_state=config['cross_validation']['random_seed'])

cv_folds = rkf.split(range(samples.shape[0]))

# convert from list of tuples to a nested dict for better snakemake compatibility
cv_folds = {'{0:02d}'.format(i + 1): { 'train': x[0], 'test': x[1] } for i, x in enumerate(cv_folds)}

# list of sequential numbers equal to the total number of folds to be tested;
# used to let snakemake know what files are to be expected
num_folds = config['cross_validation']['num_splits'] * config['cross_validation']['num_repeats']
cv_indices = [f'{x:02}' for x in list(range(1, num_folds + 1))]

#
# Rules
#
#
rule all:
    input: expand(join(output_dir, '{cv}/train/training_sets/selected/{drug}.tsv.gz'), cv=cv_indices, drug=drug_names)

#
# Model training
#
rule train_model:
    input: join(output_dir, '{cv}/train/training_sets/selected/{drug}.tsv.gz')
    output: join(output_dir, '{cv}/train/models/{drug}.rda')
    threads: config['num_threads']['train_model']
    script: 'scripts/train_model.R'

#
# Feature selection
#
if config['dimension_reduction_late']['enabled']:
    inputs = [join(output_dir, '{cv}/train/training_sets/dimension_reduced/{drug}.tsv.gz',
              join(output_dir, '{cv}/train/training_sets/dimension_reduced/{drug}_projection_matrix.rda']

else:
    inputs = [join(output_dir, '{cv}/train/training_sets/full/{drug}.tsv.gz']

rule perform_feature_selection:
    input: 
    output: join(output_dir, '{cv}/train/training_sets/selected/{drug}.tsv.gz'),
    threads: config['num_threads']['train_model']
    script:
        'scripts/select_features.R'

#
# Imputation (TODO / Optional)
#

#
# Late dimension reduction (Optional)
#
if config['dimension_reduction_late']['enabled']:
    rule reduce_training_set_dimension:
        input: join(output_dir, '{cv}/train/training_sets/full/{drug}.tsv.gz')
        output:
            join(output_dir, '{cv}/train/training_sets/dimension_reduced/{drug}.tsv.gz'),
            join(output_dir, '{cv}/train/training_sets/dimension_reduced/{drug}_projection_matrix.rda')
        script: 'scripts/reduce_dimensions_late.R'

#
# Create training set
#
rule create_training_set:
    input:
        rna=join(output_dir, '{cv}/train/features/filtered/rna.tsv.gz'),
        cnv=join(output_dir, '{cv}/train/features/filtered/cnv.tsv.gz'),
        var=join(output_dir, '{cv}/train/features/filtered/var.tsv.gz'),
        response=join(output_dir, '{cv}/train/response/{drug}.tsv.gz')
    output:
        join(output_dir, '{cv}/train/training_sets/full/{drug}.tsv.gz')
    script:
        'scripts/create_training_set.R'

#
# Feature filtering
#
if config['dimension_reduction_early']['enabled']:
    subdir = 'dimension_reduced'
else:
    if config['gene_set_projection']['enabled']:
        subdir = 'gene_set_projected'
    else:
        subdir = 'raw'

rule filter_rna_features:
    input: join(output_dir, '{{cv}}/train/features/{}/rna.tsv.gz'.format(subdir))
    output: join(output_dir, '{cv}/train/features/filtered/rna.tsv.gz')
    script: 'scripts/filter_features.R'

rule filter_cnv_features:
    input: join(output_dir, '{{cv}}/train/features/{}/cnv.tsv.gz'.format(subdir))
    output: join(output_dir, '{cv}/train/features/filtered/cnv.tsv.gz')
    script: 'scripts/filter_features.R'

rule filter_var_features:
    input: join(output_dir, '{{cv}}/train/features/{}/var.tsv.gz'.format(subdir)) 
    output: join(output_dir, '{cv}/train/features/filtered/var.tsv.gz')
    script: 'scripts/filter_features.R'

#
# Early dimension reduction (Optional)
#
if config['dimension_reduction_early']['enabled']:
    if config['gene_set_projection']['enabled']:
        subdir = 'gene_set_projected'
    else:
        subdir = 'raw'

    rule reduce_rna_dimension:
        input: join(output_dir, '{{cv}}/train/features/{}/rna.tsv.gz'.format(subdir))
        output: join(output_dir, '{cv}/train/features/dimension_reduced/rna.tsv.gz')
        script: 'scripts/reduce_dimensions_early.R'

    rule reduce_cnv_dimension:
        input: join(output_dir, '{{cv}}/train/features/{}/cnv.tsv.gz'.format(subdir))
        output: join(output_dir, '{cv}/train/features/dimension_reduced/cnv.tsv.gz')
        script: 'scripts/reduce_dimensions_early.R'

    rule reduce_var_dimension:
        input: join(output_dir, '{{cv}}/train/features/{}/var.tsv.gz'.format(subdir)) 
        output: join(output_dir, '{cv}/train/features/dimension_reduced/var.tsv.gz')
        script: 'scripts/reduce_dimensions_early.R'

#
# Gene set aggregation (Optional)
#
if config['gene_set_projection']['enabled']:
    rule project_rna_gene_sets:
        input: join(output_dir, '{cv}/train/features/raw/rna.tsv.gz')
        output: join(output_dir, '{cv}/train/features/gene_set_projected/rna.tsv.gz')
        script: 'scripts/project_gene_sets.R'

    rule project_cnv_gene_sets:
        input: join(output_dir, '{cv}/train/features/raw/cnv.tsv.gz')
        output: join(output_dir, '{cv}/train/features/gene_set_projected/cnv.tsv.gz')
        script: 'scripts/project_gene_sets.R'

    rule project_var_gene_sets:
        input: join(output_dir, '{cv}/train/features/raw/var.tsv.gz')
        output: join(output_dir, '{cv}/train/features/gene_set_projected/var.tsv.gz')
        script: 'scripts/project_gene_sets.R'

#
# Create cross validation splits
#
rule create_rna_cv_folds:
    input: join(input_dir, config['features']['rna'])
    output:
        join(output_dir, '{cv}/train/features/raw/rna.tsv.gz'),
        join(output_dir, '{cv}/test/features/raw/rna.tsv.gz'),
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_cnv_cv_folds:
    input: join(input_dir, config['features']['cnv'])
    output:
        join(output_dir, '{cv}/train/features/raw/cnv.tsv.gz'),
        join(output_dir, '{cv}/test/features/raw/cnv.tsv.gz'),
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_var_cv_folds:
    input: join(input_dir, config['features']['var'])
    output:
        join(output_dir, '{cv}/train/features/raw/var.tsv.gz'),
        join(output_dir, '{cv}/test/features/raw/var.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_response_folds:
    input: join(input_dir, 'response/{drug}.tsv.gz')
    output:
        join(output_dir, '{cv}/train/response/{drug}.tsv.gz'),
        join(output_dir, '{cv}/test/response/{drug}.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

