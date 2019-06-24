#
# Multiple Myeloma Dose Response Prediction Model Iteration
# V. Keith Hughitt
# June 2019
#
import glob
import os
from pathlib import Path

configfile: "config.yaml"

# base data directory
data_dir = os.path.join('data', config['name'])

# feature dataset versions
rna_datasets = [Path(x).name.replace('.tsv.gz', '') for x in glob.glob(os.path.join(data_dir, 'raw/features/rnaseq/*'))]
cnv_datasets = [Path(x).name.replace('.tsv.gz', '') for x in glob.glob(os.path.join(data_dir, 'raw/features/cnv/*'))]
var_datasets = [Path(x).name.replace('.tsv.gz', '') for x in glob.glob(os.path.join(data_dir, 'raw/features/variants/*'))]

# drug response input files
response_files = os.listdir(os.path.join(data_dir, 'raw/response'))
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
        [expand(os.path.join(data_dir, 'train/orig/{response}/train_RNA_{rna}_CNV_{cnv}_VAR_{var}.tsv.gz'), rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets),
         expand(os.path.join(data_dir, 'train/pca/{response}/train_RNA_{rna}_CNV_{cnv}_VAR_{var}.tsv.gz'), rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets)]
         #expand(os.path.join(data_dir, 'train/gene_set/{response}/train_{gene_set}_RNA_{rna}_CNV_{cnv}_VAR_{var}.tsv.gz'), gene_set=gene_sets, rna=rna_datasets, cnv=cnv_datasets, var=var_datasets, response=response_datasets)]

rule create_training_sets:
    input:
        os.path.join(data_dir, 'raw/features/rnaseq/{rna}.tsv.gz'),
        os.path.join(data_dir, 'raw/features/cnv/{cnv}.tsv.gz'),
        os.path.join(data_dir, 'raw/features/variants/{var}.tsv.gz'),
        os.path.join(data_dir, 'raw/response/{response}.tsv.gz')
    output:
        os.path.join(data_dir, 'train/orig/{response}/train_RNA_{rna}_CNV_{cnv}_VAR_{var}.tsv.gz')
    script:
        'scripts/create_training_set.R'

rule create_pca_projected_training_sets:
    input:
        os.path.join(data_dir, 'processed/features/rnaseq/pca/{rna}.tsv.gz'),
        os.path.join(data_dir, 'processed/features/cnv/pca/{cnv}.tsv.gz'),
        os.path.join(data_dir, 'processed/features/variants/pca/{var}.tsv.gz'),
        os.path.join(data_dir, 'raw/response/{response}.tsv.gz')
    output:
        os.path.join(data_dir, 'train/pca/{response}/train_RNA_{rna}_CNV_{cnv}_VAR_{var}.tsv.gz')
    script:
        'scripts/create_training_set.R'

rule create_gene_set_projected_training_sets:
    input:
        os.path.join(data_dir, 'processed/features/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz'),
        os.path.join(data_dir, 'processed/features/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz'),
        os.path.join(data_dir, 'processed/features/variants/gene_sets/{gene_set}/{var}.tsv.gz'),
        os.path.join(data_dir, 'raw/response/{response}.tsv.gz')
    output:
        os.path.join(data_dir, 'train/gene_set/{response}/train_{gene_set}_RNA_{rna}_CNV_{cnv}_VAR_{var}.tsv.gz')
    script:
        'scripts/create_training_set.R'

rule rna_gene_sets:
    input:
        features = os.path.join(data_dir, 'raw/features/rnaseq/{rna}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: os.path.join(data_dir, 'processed/features/rnaseq/gene_sets/{gene_set}/{rna}.tsv.gz')
    script: 'scripts/gene_set_project.R'

rule cnv_gene_sets:
    input:
        features = os.path.join(data_dir, 'raw/features/cnv/{cnv}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: os.path.join(data_dir, 'processed/features/cnv/gene_sets/{gene_set}/{cnv}.tsv.gz')
    script: 'scripts/gene_set_project.R'

rule variant_gene_sets:
    input:
        features = os.path.join(data_dir, 'raw/features/variants/{var}.tsv.gz'),
        gene_set = 'gene_sets/{gene_set}.gmt.gz'
    output: os.path.join(data_dir, 'processed/features/variants/gene_sets/{gene_set}/{var}.tsv.gz')
    script: 'scripts/gene_set_project.R'

rule rna_pca:
    input: os.path.join(data_dir, 'raw/features/rnaseq/{rna}.tsv.gz')
    output: os.path.join(data_dir, 'processed/features/rnaseq/pca/{rna}.tsv.gz')
    script: 'scripts/pca_project.R'

rule cnv_pca:
    input: os.path.join(data_dir, 'raw/features/cnv/{cnv}.tsv.gz')
    output: os.path.join(data_dir, 'processed/features/cnv/pca/{cnv}.tsv.gz')
    script: 'scripts/pca_project.R'

rule variant_pca:
    input: os.path.join(data_dir, 'raw/features/variants/{var}.tsv.gz')
    output: os.path.join(data_dir, 'processed/features/variants/pca/{var}.tsv.gz')
    script: 'scripts/pca_project.R'

