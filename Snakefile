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
data_dir = join('data', config['name'])
build_dir = join('output', config['name'], str(config['version']))

# load sample metadata
samples = pd.read_csv(join(data_dir, 'metadata', 'samples.tsv'))

# feature dataset versions
rna_datasets = [Path(x).name.replace('.tsv.gz', '') for x in glob.glob(join(data_dir, 'rnaseq/*.tsv.gz'))]
cnv_datasets = [Path(x).name.replace('.tsv.gz', '') for x in glob.glob(join(data_dir, 'cnv/*.tsv.gz'))]
var_datasets = [Path(x).name.replace('.tsv.gz', '') for x in glob.glob(join(data_dir, 'variants/*tsv.gz'))]

# drug response input files
response_files = [Path(x).name for x in glob.glob(join(data_dir, 'response/*.tsv.gz'))]
response_datasets = [x.replace('.tsv.gz', '') for x in response_files]

# gene sets
gene_set_files = [os.path.basename(x) for x in glob.glob('gene_sets/*.gmt.gz')]
gene_sets = [Path(x).name.replace('.gmt.gz', '') for x in gene_set_files]

#
# Cross-validation setup
#
rkf = RepeatedKFold(n_splits=config['cv']['num_splits'],
                    n_repeats=config['cv']['num_repeats'],
                    random_state=config['cv']['random_seed'])

cv_folds = rkf.split(range(samples.shape[0]))

# convert from list of tuples to a nested dict for better snakemake compatibility
cv_folds = {i: { 'train': x[0], 'test': x[1] } for i, x in enumerate(cv_folds)}

# list of sequential numbers equal to the total number of folds to be tested;
# used to let snakemake know what files are to be expected
cv_indices = list(range(config['cv']['num_splits'] * config['cv']['num_repeats']))

#
# Rules
#
#
rule all:
    input:
        train=[expand(join(build_dir, '{cv}/train/combined/orig/{response}/{rna}/{cnv}/{var}/combined.tsv.gz'),
                      cv=cv_indices, rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets)],
               #expand(join(build_dir, '{cv}/train/combined/pca/{response}/{rna}/{cnv}/{var}/combined.tsv.gz'),
               #       cv=cv_indices, rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets),
               #expand(join(build_dir, '{cv}/train/combined/gene_set_{gene_set}/{response}/{rna}/{cnv}/{var}/combined.tsv.gz'),
               #       cv=cv_indices, gene_set=gene_sets, rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets)],
        test=[expand(join(build_dir, '{cv}/test/raw/rnaseq/{rna}.tsv.gz'), cv=cv_indices, rna=rna_datasets),
              expand(join(build_dir, '{cv}/test/raw/cnv/{cnv}.tsv.gz'), cv=cv_indices, cnv=cnv_datasets),
              expand(join(build_dir, '{cv}/test/raw/variants/{var}.tsv.gz'), cv=cv_indices, var=var_datasets),
              expand(join(build_dir, '{cv}/test/raw/response/{response}.tsv.gz'), cv=cv_indices, response=response_datasets)]

#
# combine feature and response data
#
rule create_training_sets:
    input:
        join(build_dir, '{cv}/train/filtered/rnaseq/orig/{rna}.tsv.gz'),
        join(build_dir, '{cv}/train/filtered/cnv/orig/{cnv}.tsv.gz'),
        join(build_dir, '{cv}/train/filtered/variants/orig/{var}.tsv.gz'),
        join(build_dir, '{cv}/train/raw/response/{response}.tsv.gz')
    output:
        join(build_dir, '{cv}/train/combined/orig/{response}/{rna}/{cnv}/{var}/train.tsv.gz')
    script:
        'scripts/create_training_set.R'

rule create_pca_projected_training_sets:
    input:
        join(build_dir, '{cv}/train/filtered/rnaseq/pca/{rna}.tsv.gz'),
        join(build_dir, '{cv}/train/filtered/cnv/pca/{cnv}.tsv.gz'),
        join(build_dir, '{cv}/train/filtered/variants/pca/{var}.tsv.gz'),
        join(build_dir, '{cv}/train/raw/response/{response}.tsv.gz')
    output:
        join(build_dir, '{cv}/train/combined/pca/{response}/{rna}/{cnv}/{var}/train.tsv.gz')
    script:
        'scripts/create_training_set.R'

rule create_gene_set_projected_training_sets:
    input:
        join(build_dir, '{cv}/train/filtered/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz'),
        join(build_dir, '{cv}/train/filtered/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz'),
        join(build_dir, '{cv}/train/filtered/variants/gene_sets/{gene_set}/{var}.tsv.gz'),
        join(build_dir, '{cv}/train/raw/response/{response}.tsv.gz')
    output:
        join(build_dir, '{cv}/train/combined/gene_set_{gene_set}/{response}/{rna}/{cnv}/{var}/train.tsv.gz')
    script:
        'scripts/create_training_set.R'

#
# Feature selection
#
rule select_rna_features:
    input: join(build_dir, '{cv}/train/raw/rnaseq/{rna}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/rnaseq/orig/{rna}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_cnv_features:
    input: join(build_dir, '{cv}/train/raw/cnv/{cnv}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/cnv/orig/{cnv}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_variant_features:
    input: join(build_dir, '{cv}/train/raw/variants/{var}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/variants/orig/{var}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_rna_pca_features:
    input: join(build_dir, '{cv}/train/processed/rnaseq/pca/{rna}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/rnaseq/pca/{rna}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_cnv_pca_features:
    input: join(build_dir, '{cv}/train/processed/cnv/pca/{cnv}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/cnv/pca/{cnv}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_variant_pca_features:
    input: join(build_dir, '{cv}/train/processed/variants/pca/{var}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/variants/pca/{var}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_rna_gene_set_features:
    input: join(build_dir, '{cv}/train/processed/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_cnv_gene_set_features:
    input: join(build_dir, '{cv}/train/processed/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_variant_gene_set_features:
    input: join(build_dir, '{cv}/train/processed/variants/gene_sets/{gene_set}/{var}.tsv.gz')
    output: join(build_dir, '{cv}/train/filtered/variants/gene_sets/{gene_set}/{var}.tsv.gz')
    script: 'scripts/select_features.R'

#
# Gene set aggregation
#
rule rna_gene_sets:
    input:
        features = join(build_dir, '{cv}/train/raw/rnaseq/{rna}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: join(build_dir, '{cv}/train/processed/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz')
    script: 'scripts/project_gene_sets.R'

rule cnv_gene_sets:
    input:
        features = join(build_dir, '{cv}/train/raw/cnv/{cnv}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: join(build_dir, '{cv}/train/processed/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz')
    script: 'scripts/project_gene_sets.R'

rule variant_gene_sets:
    input:
        features = join(build_dir, '{cv}/train/raw/variants/{var}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: join(build_dir, '{cv}/train/processed/variants/gene_sets/{gene_set}/{var}.tsv.gz')
    script: 'scripts/project_gene_sets.R'

#
# PCA projection
#
rule rna_pca:
    input: join(build_dir, '{cv}/train/raw/rnaseq/{rna}.tsv.gz')
    output: join(build_dir, '{cv}/train/processed/rnaseq/pca/{rna}.tsv.gz')
    script: 'scripts/project_pca.R'

rule cnv_pca:
    input: join(build_dir, '{cv}/train/raw/cnv/{cnv}.tsv.gz')
    output: join(build_dir, '{cv}/train/processed/cnv/pca/{cnv}.tsv.gz')
    script: 'scripts/project_pca.R'

rule variant_pca:
    input: join(build_dir, '{cv}/train/raw/variants/{var}.tsv.gz')
    output: join(build_dir, '{cv}/train/processed/variants/pca/{var}.tsv.gz')
    script: 'scripts/project_pca.R'

#
# Cross validation splits
#
rule create_rna_cv_folds:
    input: join(data_dir, 'rnaseq/{rna}.tsv.gz')
    output:
        join(build_dir, '{cv}/train/raw/rnaseq/{rna}.tsv.gz'),
        join(build_dir, '{cv}/test/raw/rnaseq/{rna}.tsv.gz'),
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_cnv_cv_folds:
    input: join(data_dir, 'cnv/{cnv}.tsv.gz')
    output:
        join(build_dir, '{cv}/train/train/raw/cnv/{cnv}.tsv.gz'),
        join(build_dir, '{cv}/test/train/raw/cnv/{cnv}.tsv.gz'),
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_variant_folds:
    input: join(data_dir, 'variants/{var}.tsv.gz')
    output:
        join(build_dir, '{cv}/train/train/raw/variants/{var}.tsv.gz'),
        join(build_dir, '{cv}/test/train/raw/variants/{var}.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_response_folds:
    input: join(data_dir, 'response/{response}.tsv.gz')
    output:
        join(build_dir, '{cv}/train/raw/response/{response}.tsv.gz'),
        join(build_dir, '{cv}/test/raw/response/{response}.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_response_cv_folds.R'

