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
import numpy as np
import pandas as pd
from os.path import join
from pathlib import Path
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.cluster import KMeans

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
if config['cross_validation']['num_splits'] == 1:
    #
    # No cross-validation
    #
    cv_folds = ["alldata"]
    cv_indices = ["alldata"]
    model_save = "rds"
else:
    #
    # Stratified Repeated cross validation
    #
    model_save = "tsv.gz"

    rskf = RepeatedStratifiedKFold(n_splits=config['cross_validation']['num_splits'],
                                   n_repeats=config['cross_validation']['num_repeats'],
                                   random_state=config['cross_validation']['random_seed'])

    # create a response-index dict to store CV indices
    cv_folds = {}

    # for each drug, cluster response data into two groups, if possible, and use to
    # generate balanced (stratefied) CV splits
    for response_file in response_files:
        # load response data and convert to an n x 1 array
        response = pd.read_csv(join(input_dir, 'response', response_file), sep='\t')
        arr = response.iloc[0, 1:].to_numpy()
        dat = arr.reshape(-1, 1)

        # perform kmeans clustering
        clusters = KMeans(n_clusters=2, random_state=0).fit(dat).labels_

        # if too few clusters are found, use extremes instead
        group_sizes = np.bincount(clusters)

        # group_sizes     
        # array([29, 14])

        min_group_size = config['cross_validation']['min_group_size']

        if min(group_sizes) < min_group_size:
            # determine which group is associated with lower values
            if np.mean(dat[clusters == 0]) < np.mean(dat[clusters == 1]):
                small_value_group = 0
            else:
                small_value_group = 1
            
            # determine which group is smaller, and extend it to the minimum size
            if ((group_sizes[0] <= group_sizes[1] and small_value_group == 0) or 
                (group_sizes[1] <= group_sizes[0] and small_value_group == 1)):
                # if smaller group contains the lower values, extend to the N smallest values
                index = arr.argsort()[:min_group_size]
            else:
                # otherwise, if smaller group contains the higher values, extend to
                # the N largest values
                index = arr.argsort()[-min_group_size:]

            clusters = np.zeros(len(clusters), dtype=np.int8)
            clusters[index] = 1

        # store CV fold indices
        drug_folds = rskf.split(range(samples.shape[0]), clusters)

        # convert from list of tuples to a nested dict for better snakemake/R compatibility
        drug_folds = {'{0:02d}'.format(i + 1): { 'train': x[0], 'test': x[1] } for i, x in enumerate(drug_folds)}

        # list of sequential numbers equal to the total number of folds to be tested;
        # used to let snakemake know what files are to be expected
        num_folds = config['cross_validation']['num_splits'] * config['cross_validation']['num_repeats']
        cv_indices = [f'{x:02}' for x in list(range(1, num_folds + 1))]

        cv_folds[response_file.replace('.tsv.gz', '')] = drug_folds

# specify which rules are run locally
localrules: all, 
    create_training_set, 
    create_rna_cv_folds, 
    create_cnv_cv_folds,
    create_var_cv_folds,
    create_response_folds

data_aggregations = config['model_combinations']['aggregations']
models = config['model_combinations']['models']
dim_reducts = config['model_combinations']['dim_reducts']
feat_select = config['model_combinations']['feat_select']

# specify format of wildcards to prevent ambiguous file names
wildcard_constraints:
   cv = "(\d+)|(alldata)",
   drug = "[^\/]+", 
   aggregation = "({})".format(")|(".join(data_aggregations)), 
   feat = "({})".format(")|(".join(feat_select)), 
   dimreduct = "({})".format(")|(".join(dim_reducts)), 
   model = "({})".format(")|(".join(models))


#
# Rules
#
#
rule all:
    input: expand(join(output_dir, '{{aggregation}}/{{drug}}/{{cv}}/models/{{feat}}/{{dimreduct}}/{{model}}.{}'.format(model_save)), aggregation=data_aggregations, dimreduct=dim_reducts, model=models, cv=cv_indices, drug=drug_names, feat=feat_select)

#
# Model training
#

rule evaluate_model:
    input: join(output_dir, '{aggregation}/{drug}/{cv}/training_sets/dimension_reduced/{dimreduct}/{feat}/response.tsv.gz')
    output: join(output_dir, '{{aggregation}}/{{drug}}/{{cv}}/models/{{feat}}/{{dimreduct}}/{{model}}.{}'.format(model_save))
    threads: config['num_threads']['train_model']
    script: 'scripts/eval_model.R'

#
# Late dimension reduction 
#
rule reduce_training_set_dimension:
    input: join(output_dir, '{aggregation}/{drug}/{cv}/training_sets/selected/{feat}/response.tsv.gz')
    output: 
        join(output_dir, '{aggregation}/{drug}/{cv}/training_sets/dimension_reduced/{dimreduct}/{feat}/response.tsv.gz'),
        join(output_dir, '{aggregation}/{drug}/{cv}/training_sets/dimension_reduced/{dimreduct}/{feat}/extra.tsv.gz') 
    script: 'scripts/reduce_dimensions_late.R'

#
# Feature selection
#
rule perform_feature_selection:
    input: join(output_dir, '{aggregation}/{drug}/{cv}/training_sets/full/response.tsv.gz')
    output: join(output_dir, '{aggregation}/{drug}/{cv}/training_sets/selected/{feat}/response.tsv.gz'),
    threads: config['num_threads']['feature_selection']
    script:
        'scripts/select_features.R'
 

#
# Create training set
#
rule create_training_set:
    input:
        rna=join(output_dir, '{aggregation}/{drug}/{cv}/features/filtered/rna.tsv.gz'),
        cnv=join(output_dir, '{aggregation}/{drug}/{cv}/features/filtered/cnv.tsv.gz'),
        var=join(output_dir, '{aggregation}/{drug}/{cv}/features/filtered/var.tsv.gz'),
        response=join(output_dir, '{aggregation}/{drug}/{cv}/response/response.tsv.gz')
    output:
        join(output_dir, '{aggregation}/{drug}/{cv}/training_sets/full/response.tsv.gz')
    script:
        'scripts/create_training_set.R'


#
# Feature filtering
#
rule filter_rna_features:
    input: join(output_dir, '{aggregation}/{drug}/{cv}/features/rna.tsv.gz')
    output: join(output_dir, '{aggregation}/{drug}/{cv}/features/filtered/rna.tsv.gz')
    script: 'scripts/filter_features.R'

rule filter_cnv_features:
    input: join(output_dir, '{aggregation}/{drug}/{cv}/features/cnv.tsv.gz')
    output: join(output_dir, '{aggregation}/{drug}/{cv}/features/filtered/cnv.tsv.gz')
    script: 'scripts/filter_features.R'

rule filter_var_features:
    input: join(output_dir, '{aggregation}/{drug}/{cv}/features/var.tsv.gz') 
    output: join(output_dir, '{aggregation}/{drug}/{cv}/features/filtered/var.tsv.gz')
    script: 'scripts/filter_features.R'

#
# Create cross validation splits
#

rule create_rna_cv_folds:
    input: join(output_dir, '{aggregation}/{drug}/features/rna.tsv.gz')
    output: join(output_dir, '{aggregation}/{drug}/{cv}/features/rna.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_cnv_cv_folds:
    input: join(output_dir, '{aggregation}/{drug}/features/cnv.tsv.gz')
    output: join(output_dir, '{aggregation}/{drug}/{cv}/features/cnv.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_var_cv_folds:
    input: join(output_dir, '{aggregation}/{drug}/features/var.tsv.gz')
    output: join(output_dir, '{aggregation}/{drug}/{cv}/features/var.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

rule create_response_folds:
    input: join(input_dir, 'response/{drug}.tsv.gz')
    output: join(output_dir, '{aggregation}/{drug}/{cv}/response/response.tsv.gz')
    params:
        cv_folds=cv_folds
    script: 'scripts/create_cv_folds.R'

#
# Feature aggregation
#
rule aggregate_rna:
    input: join(input_dir, config['features']['rna'])
    output: join(output_dir, '{aggregation}/{drug}/features/rna.tsv.gz')
    script: 'scripts/aggregate_features.R'

rule aggregate_cnv:
    input: join(input_dir, config['features']['cnv'])
    output: join(output_dir, '{aggregation}/{drug}/features/cnv.tsv.gz')
    script: 'scripts/aggregate_features.R'

rule aggregate_var:
    input: join(input_dir, config['features']['var'])
    output: join(output_dir, '{aggregation}/{drug}/features/var.tsv.gz')
    script: 'scripts/aggregate_features.R'


