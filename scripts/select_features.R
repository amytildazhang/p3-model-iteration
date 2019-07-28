#!/bin/env Rscript
#
# Feature selection
#
library(tidyverse)
source('scripts/feature_selection_methods.R')

set.seed(1)

# load full training set
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# drop sample ids and convert to a data frame
dat <- as.data.frame(dat)
sample_ids <- dat[, 1]

# drop samples with missing response values
dat <- dat[!is.na(dat$response), ]

# drop samples with response values of "Inf" (e.g. GDSC)
mask <- dat$response != Inf
dat <- dat[mask, ]

# drop any features with missing values
dat <- dat[, complete.cases(t(dat))]

# if the number of features in input dataset is already less than or equal to the
# minimum desired number of features, simply create a copy of the data and exit
if (nrow(dat) <= snakemake@config$feature_selection$min_features) {
  message("Number of features remaining is already at or below desired number; skipping feature selection...")
  write_tsv(dat, snakemake@output[[1]])
  quit(save = 'no')
}

# non-feature column indices
CV_IND <- 1
SAMPLE_IND <- 2 

# separate out CV indicator
train_indices <- dat[, CV_IND] == 1

METHOD <- snakemake@wildcards$feat

# perform feature selection (first attempt)
if (METHOD == 'boruta') {
  features <- boruta_feature_selection(dat[train_indices, -c(CV_IND, SAMPLE_IND)], snakemake) 
} else if (METHOD == 'rfe') {
  features <- rfe_feature_selection(dat[train_indices, -c(CV_IND, SAMPLE_IND)], snakemake) 
} else if (METHOD == 'distance') {
  features <- dcor_feature_selection(dat[train_indices, -c(CV_IND, SAMPLE_IND)], snakemake)
}  else if (METHOD == 'none') {
  # if the feature selection method is set to "none", we can stop here and
  # simply return the full dataset
  write_tsv(dat, snakemake@output[[1]])
  quit(save = 'no')
} else {
  stop("Invalid feature selection method specified!")
}

message(sprintf("Found %d features during first round of selection using %s...",
                length(features),
                snakemake@config$feature_selection$method))
print(features)

# if too few features found, attempt fallback method
if (length(features) < snakemake@config$feature_selection$min_features) {
  message(sprintf("Insufficient features found using %s! Falling back on %s...",
                  snakemake@config$feature_selection$method,
                  snakemake@config$feature_selection$fallback))
  # fallback: boruta
  if (snakemake@config$feature_selection$fallback == 'boruta') {
    features <- boruta_feature_selection(dat[train_indices, -c(CV_IND, SAMPLE_IND)], snakemake) 
  } else if (snakemake@config$feature_selection$fallback == 'rfe') {
    # fallback: rfe
    features <- rfe_feature_selection(dat[train_indices, -c(CV_IND, SAMPLE_IND)], snakemake) 
  } else if (snakemake@config$feature_selection$fallback == 'distance') {
    # fallback: rfe
    features <- dcor_feature_selection(dat[train_indices, -c(CV_IND, SAMPLE_IND)], snakemake)
  } else if (snakemake@config$feature_selection$fallback == 'none') {
    # if the feature selection method is set to "none", we can stop here and
    # simply return the full dataset
    write_tsv(dat, snakemake@output[[1]])
    quit(save = 'no')
  } else {
    stop("Invalid fallback feature selection method specified!")
  }

  message(sprintf("Found %d features during fallback round of selection using %s...",
                  length(features),
                  snakemake@config$feature_selection$fallback))
}

# remove unselected features
dat <- dat[, colnames(dat) %in% c('sample_id', 'symbol', 'train_indices', features, 'response')]

# store result
write_tsv(dat, snakemake@output[[1]])


