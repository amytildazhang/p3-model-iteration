#
# Dose Response Prediction Model Iteration
# V. Keith Hughitt, Amy Zhang
# June 2019
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
# train models
#
#rule train_rf_models:
#    input: join(output_dir, '{cv}/train/training_sets/orig/{response}.tsv.gz')
#    output: join(output_dir, '{cv}/train/models/orig/{response}.tsv.gz')
#    script: 'scripts/train_random_forest_model.R'


#
# combine feature and response data
#
#rule train_rf_models:
#    input: join(output_dir, '{cv}/train/training_sets/orig/{response}.tsv.gz')
#    output: join(output_dir, '{cv}/train/models/orig/{response}.tsv.gz')
#    script: 'scripts/train_random_forest_model.R'


#
# Feature selection
#
rule select_rna_features:
    input: join(output_dir, '{cv}/train/raw/orig/rna.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/orig/rna.tsv.gz')
    script: 'scripts/filter_features.R'

rule select_cnv_features:
    input: join(output_dir, '{cv}/train/raw/orig/cnv.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/orig/cnv.tsv.gz')
    script: 'scripts/filter_features.R'

rule select_var_features:
    input: join(output_dir, '{cv}/train/raw/orig/var.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/orig/var.tsv.gz')
    script: 'scripts/filter_features.R'

#rule select_rna_pca_features:
#    input: join(output_dir, '{cv}/train/raw/pca_projected/rna.tsv.gz')
#    output: join(output_dir, '{cv}/train/filtered/pca_projected/rna.tsv.gz')
#    script: 'scripts/filter_features.R'

#rule select_cnv_pca_features:
#    input: join(output_dir, '{cv}/train/raw/pca_projected/cnv/pca/{cnv}.tsv.gz')
#    output: join(output_dir, '{cv}/train/filtered/pca_projected/cnv/pca/{cnv}.tsv.gz')
#    script: 'scripts/filter_features.R'

#rule select_var_pca_features:
#    input: join(output_dir, '{cv}/train/raw/pca_projected/var.tsv.gz')
#    output: join(output_dir, '{cv}/train/filtered/pca_projected/var.tsv.gz')
#    script: 'scripts/filter_features.R'

rule select_rna_gene_set_features:
    input: join(output_dir, '{cv}/train/raw/gene_set_projected/rna.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/gene_set_projected/rna.tsv.gz')
    script: 'scripts/filter_features.R'

rule select_cnv_gene_set_features:
    input: join(output_dir, '{cv}/train/raw/gene_set_projected/cnv.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/gene_set_projected/cnv.tsv.gz')
    script: 'scripts/filter_features.R'

rule select_var_gene_set_features:
    input: join(output_dir, '{cv}/train/raw/gene_set_projected/var.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/gene_set_projected/var.tsv.gz')
    script: 'scripts/filter_features.R'

#
# Model training (TODO)
#
#rule train_rf_models:
#    input: join(output_dir, '{cv}/train/training_sets/selected/{drug}.tsv.gz')
#    output: join(output_dir, '{cv}/train/models/{drug}.tsv.gz')
#    script: 'scripts/train_random_forest_model.R'

rule perform_feature_selection:
    input: join(output_dir, '{cv}/train/training_sets/full/{drug}.tsv.gz')
    output: join(output_dir, '{cv}/train/training_sets/selected/{drug}.tsv.gz')
    script:
        'scripts/select_features.R'

#
# Imputation (TODO / Optional)
#

#
# Create training set
#
rule create_training_sets:
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
if config['pca_projection']['enabled']:
    subdir = 'pca_projected'
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
# PCA projection (Optional)
#
if config['pca_projection']['enabled']:
    if config['gene_set_projection']['enabled']:
        subdir = 'gene_set_projected'
    else:
        subdir = 'raw'

    rule project_rna_pca:
        input: join(output_dir, '{{cv}}/train/features/{}/rna.tsv.gz'.format(subdir))
        output: join(output_dir, '{cv}/train/features/pca_projected/rna.tsv.gz')
        script: 'scripts/project_pca.R'
    rule project_cnv_pca:
        input: join(output_dir, '{{cv}}/train/features/{}/cnv.tsv.gz'.format(subdir))
        output: join(output_dir, '{cv}/train/features/pca_projected/cnv.tsv.gz')
        script: 'scripts/project_pca.R'
    rule project_var_pca:
        input: join(output_dir, '{{cv}}/train/features/{}/var.tsv.gz'.format(subdir)) 
        output: join(output_dir, '{cv}/train/features/pca_projected/var.tsv.gz')
        script: 'scripts/project_pca.R'

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

