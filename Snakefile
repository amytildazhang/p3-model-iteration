#
# Multiple Myeloma Dose Response Prediction Model Iteration
# V. Keith Hughitt
# June 2019
#
import glob
import os
from os.path import join
from pathlib import Path

configfile: "config.yaml"

# base data directory
data_dir = join('data', config['name'])
build_dir = join('output', config['name'], str(config['version']))

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
# Rules
#
#
rule all:
    input:
        [expand(join(build_dir, 'train/orig/{response}/{rna}/{cnv}/{var}/train.tsv.gz'), rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets),
         expand(join(build_dir, 'train/pca/{response}/{rna}/{cnv}/{var}/train.tsv.gz'),  rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets),
         expand(join(build_dir, 'train/gene_set_{gene_set}/{response}/{rna}/{cnv}/{var}/train.tsv.gz'), gene_set=gene_sets, rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets)]

#
# Training set construction
#
rule create_training_sets:
    input:
        join(build_dir, 'filtered/rnaseq/orig/{rna}.tsv.gz'),
        join(build_dir, 'filtered/cnv/orig/{cnv}.tsv.gz'),
        join(build_dir, 'filtered/variants/orig/{var}.tsv.gz'),
        join(data_dir, 'response/{response}.tsv.gz')
    output:
        join(build_dir, 'train/orig/{response}/{rna}/{cnv}/{var}/train.tsv.gz')
    script:
        'scripts/create_training_set.R'

rule create_pca_projected_training_sets:
    input:
        join(build_dir, 'filtered/rnaseq/pca/{rna}.tsv.gz'),
        join(build_dir, 'filtered/cnv/pca/{cnv}.tsv.gz'),
        join(build_dir, 'filtered/variants/pca/{var}.tsv.gz'),
        join(data_dir, 'response/{response}.tsv.gz')
    output:
        join(build_dir, 'train/pca/{response}/{rna}/{cnv}/{var}/train.tsv.gz')
    script:
        'scripts/create_training_set.R'

rule create_gene_set_projected_training_sets:
    input:
        join(build_dir, 'filtered/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz'),
        join(build_dir, 'filtered/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz'),
        join(build_dir, 'filtered/variants/gene_sets/{gene_set}/{var}.tsv.gz'),
        join(data_dir, 'response/{response}.tsv.gz')
    output:
        join(build_dir, 'train/gene_set_{gene_set}/{response}/{rna}/{cnv}/{var}/train.tsv.gz')
    script:
        'scripts/create_training_set.R'

#
# Feature selection
#
rule select_rna_features:
    input: join(data_dir, 'rnaseq/{rna}.tsv.gz')
    output: join(build_dir, 'filtered/rnaseq/orig/{rna}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_cnv_features:
    input: join(data_dir, 'cnv/{cnv}.tsv.gz')
    output: join(build_dir, 'filtered/cnv/orig/{cnv}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_variant_features:
    input: join(data_dir, 'variants/{var}.tsv.gz')
    output: join(build_dir, 'filtered/variants/orig/{var}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_rna_pca_features:
    input: join(build_dir, 'processed/rnaseq/pca/{rna}.tsv.gz')
    output: join(build_dir, 'filtered/rnaseq/pca/{rna}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_cnv_pca_features:
    input: join(build_dir, 'processed/cnv/pca/{cnv}.tsv.gz')
    output: join(build_dir, 'filtered/cnv/pca/{cnv}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_variant_pca_features:
    input: join(build_dir, 'processed/variants/pca/{var}.tsv.gz')
    output: join(build_dir, 'filtered/variants/pca/{var}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_rna_gene_set_features:
    input: join(build_dir, 'processed/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz')
    output: join(build_dir, 'filtered/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_cnv_gene_set_features:
    input: join(build_dir, 'processed/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz')
    output: join(build_dir, 'filtered/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz')
    script: 'scripts/select_features.R'

rule select_variant_gene_set_features:
    input: join(build_dir, 'processed/variants/gene_sets/{gene_set}/{var}.tsv.gz')
    output: join(build_dir, 'filtered/variants/gene_sets/{gene_set}/{var}.tsv.gz')
    script: 'scripts/select_features.R'

#
# Gene set aggregation
#
rule rna_gene_sets:
    input:
        features = join(data_dir, 'rnaseq/{rna}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: join(build_dir, 'processed/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz')
    script: 'scripts/project_gene_sets.R'

rule cnv_gene_sets:
    input:
        features = join(data_dir, 'cnv/{cnv}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: join(build_dir, 'processed/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz')
    script: 'scripts/project_gene_sets.R'

rule variant_gene_sets:
    input:
        features = join(data_dir, 'variants/{var}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: join(build_dir, 'processed/variants/gene_sets/{gene_set}/{var}.tsv.gz')
    script: 'scripts/project_gene_sets.R'

#
# PCA projection
#
rule rna_pca:
    input: join(data_dir, 'rnaseq/{rna}.tsv.gz')
    output: join(build_dir, 'processed/rnaseq/pca/{rna}.tsv.gz')
    script: 'scripts/project_pca.R'

rule cnv_pca:
    input: join(data_dir, 'cnv/{cnv}.tsv.gz')
    output: join(build_dir, 'processed/cnv/pca/{cnv}.tsv.gz')
    script: 'scripts/project_pca.R'

rule variant_pca:
    input: join(data_dir, 'variants/{var}.tsv.gz')
    output: join(build_dir, 'processed/variants/pca/{var}.tsv.gz')
    script: 'scripts/project_pca.R'

