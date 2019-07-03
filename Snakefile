#
# Multiple Myeloma Dose Response Prediction Model Iteration
# V. Keith Hughitt
# June 2019
#
# Note 2019/06/25: possible simplifications
#
#  - assume one version of each of the feature datasets
#  - create a single gene set proj. version of each data using all gene sets
#
import glob
import os
import pandas as pd
from os.path import join
from pathlib import Path
from sklearn.model_selection import RepeatedKFold

configfile: "config.yaml"

# base data directory
input_dir = join('data', config['name'])
output_dir = join('output', config['name'], str(config['version']))

# load sample metadata
samples = pd.read_csv(join(input_dir, 'metadata', 'samples.tsv'))

# gene sets
gene_set_files = [os.path.basename(x) for x in glob.glob('gene_sets/*.gmt.gz')]
gene_sets = [Path(x).name.replace('.gmt.gz', '') for x in gene_set_files]

# drug response input files
response_files = [Path(x).name for x in glob.glob(join(input_dir, 'response/*.tsv.gz'))]
drug_names = [x.replace('.tsv.gz', '') for x in response_files] 

#
# Cross-validation setup
#
rkf = RepeatedKFold(n_splits=config['cv']['num_splits'],
                    n_repeats=config['cv']['num_repeats'],
                    random_state=config['cv']['random_seed'])

cv_folds = rkf.split(range(samples.shape[0]))

# convert from list of tuples to a nested dict for better snakemake compatibility
cv_folds = {'{0:02d}'.format(i + 1): { 'train': x[0], 'test': x[1] } for i, x in enumerate(cv_folds)}

# list of sequential numbers equal to the total number of folds to be tested;
# used to let snakemake know what files are to be expected
num_folds = config['cv']['num_splits'] * config['cv']['num_repeats']
cv_indices = [f'{x:02}' for x in list(range(1, num_folds + 1))]

#
# Rules
#
#
rule all:
    input:
        expand(join(output_dir, '{cv}/train/combined/orig/{response}/train.tsv.gz'), cv=cv_indices, response=drug_names),
        expand(join(output_dir, '{cv}/train/combined/gene_set_projected/{response}/train.tsv.gz'), cv=cv_indices, response=drug_names)

#
# combine feature and response data
#
rule create_training_sets:
    input:
        join(output_dir, '{cv}/train/filtered/orig/rna.tsv.gz'),
        join(output_dir, '{cv}/train/filtered/orig/cnv.tsv.gz'),
        join(output_dir, '{cv}/train/filtered/orig/var.tsv.gz'),
        join(output_dir, '{cv}/train/raw/response/{response}.tsv.gz')
    output:
        join(output_dir, '{cv}/train/combined/orig/{response}/train.tsv.gz')
    script:
        'scripts/create_training_set.R'

#rule create_pca_projected_training_sets:
#    input:
#        join(output_dir, '{cv}/train/filtered/rna.tsv.gz'),
#        join(output_dir, '{cv}/train/filtered/cnv.tsv.gz'),
#        join(output_dir, '{cv}/train/filtered/var.tsv.gz'),
#        join(output_dir, '{cv}/train/raw/response/{response}.tsv.gz')
#    output:
#        join(output_dir, '{cv}/train/combined/pca/{response}/{rna}/{cnv}/{var}/train.tsv.gz')
#    script:
#        'scripts/create_training_set.R'

rule create_gene_set_projected_training_sets:
    input:
        join(output_dir, '{cv}/train/filtered/gene_set_projected/rna.tsv.gz'),
        join(output_dir, '{cv}/train/filtered/gene_set_projected/cnv.tsv.gz'),
        join(output_dir, '{cv}/train/filtered/gene_set_projected/var.tsv.gz'),
        join(output_dir, '{cv}/train/raw/response/{response}.tsv.gz')
    output:
        join(output_dir, '{cv}/train/combined/gene_set_projected/{response}/train.tsv.gz')
    script:
        'scripts/create_training_set.R'

#
# Feature selection
#
rule select_rna_features:
    input: join(output_dir, '{cv}/train/raw/rna.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/orig/rna.tsv.gz')
    script: 'scripts/select_features.R'

rule select_cnv_features:
    input: join(output_dir, '{cv}/train/raw/cnv.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/orig/cnv.tsv.gz')
    script: 'scripts/select_features.R'

rule select_var_features:
    input: join(output_dir, '{cv}/train/raw/var.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/orig/var.tsv.gz')
    script: 'scripts/select_features.R'

#rule select_rna_pca_features:
#    input: join(output_dir, '{cv}/train/pca_projected/rna.tsv.gz')
#    output: join(output_dir, '{cv}/train/filtered/rna.tsv.gz')
#    script: 'scripts/select_features.R'

#rule select_cnv_pca_features:
#    input: join(output_dir, '{cv}/train/pca_projected/cnv/pca/{cnv}.tsv.gz')
#    output: join(output_dir, '{cv}/train/filtered/cnv/pca/{cnv}.tsv.gz')
#    script: 'scripts/select_features.R'

#rule select_var_pca_features:
#    input: join(output_dir, '{cv}/train/pca_projected/var.tsv.gz')
#    output: join(output_dir, '{cv}/train/filtered/var.tsv.gz')
#    script: 'scripts/select_features.R'

rule select_rna_gene_set_features:
    input: join(output_dir, '{cv}/train/gene_set_projected/rna.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/gene_set_projected/rna.tsv.gz')
    script: 'scripts/select_features.R'

rule select_cnv_gene_set_features:
    input: join(output_dir, '{cv}/train/gene_set_projected/cnv.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/gene_set_projected/cnv.tsv.gz')
    script: 'scripts/select_features.R'

rule select_var_gene_set_features:
    input: join(output_dir, '{cv}/train/gene_set_projected/var.tsv.gz')
    output: join(output_dir, '{cv}/train/filtered/gene_set_projected/var.tsv.gz')
    script: 'scripts/select_features.R'

#
# Gene set aggregation
#
rule project_rna_gene_sets:
    input: join(output_dir, '{cv}/train/raw/rna.tsv.gz'),
    output: join(output_dir, '{cv}/train/gene_set_projected/rna.tsv.gz')
    script: 'scripts/project_gene_sets.R'

rule project_cnv_gene_sets:
    input: join(output_dir, '{cv}/train/raw/cnv.tsv.gz'),
    output: join(output_dir, '{cv}/train/gene_set_projected/cnv.tsv.gz')
    script: 'scripts/project_gene_sets.R'

rule project_var_gene_sets:
    input: join(output_dir, '{cv}/train/raw/var.tsv.gz'),
    output: join(output_dir, '{cv}/train/gene_set_projected/var.tsv.gz')
    script: 'scripts/project_gene_sets.R'

#
# PCA projection
#
#rule project_rna_pca:
#    input: join(output_dir, '{cv}/train/raw/rna.tsv.gz')
#    output: join(output_dir, '{cv}/train/pca_projected/rna.tsv.gz')
#    script: 'scripts/project_pca.R'

#rule project_cnv_pca:
#    input: join(output_dir, '{cv}/train/raw/cnv.tsv.gz')
#    output: join(output_dir, '{cv}/train/pca_projected/cnv.tsv.gz')
#    script: 'scripts/project_pca.R'

#rule project_var_pca:
#    input: join(output_dir, '{cv}/train/raw/var.tsv.gz')
#    output: join(output_dir, '{cv}/train/pca_projected/var.tsv.gz')
#    script: 'scripts/project_pca.R'

#
# Cross validation splits
#
rule create_rna_cv_folds:
    input: join(input_dir, config['features']['rna'])
    output:
        join(output_dir, '{cv}/train/raw/rna.tsv.gz'),
        join(output_dir, '{cv}/test/raw/rna.tsv.gz'),
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_cnv_cv_folds:
    input: join(input_dir, config['features']['cnv'])
    output:
        join(output_dir, '{cv}/train/raw/cnv.tsv.gz'),
        join(output_dir, '{cv}/test/raw/cnv.tsv.gz'),
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_var_cv_folds:
    input: join(input_dir, config['features']['var'])
    output:
        join(output_dir, '{cv}/train/raw/var.tsv.gz'),
        join(output_dir, '{cv}/test/raw/var.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_response_folds:
    input: join(input_dir, 'response/{response}.tsv.gz')
    output:
        join(output_dir, '{cv}/train/raw/response/{response}.tsv.gz'),
        join(output_dir, '{cv}/test/raw/response/{response}.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

