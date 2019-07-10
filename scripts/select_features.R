#!/bin/env Rscript
#
# Feature selection
#
library(readr)
source('scripts/feature_selection_methods.R')

set.seed(1)

# load full training set
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# drop sample ids and convert to a data frame
dat <- as.data.frame(dat)
sample_ids <- dat[, 1]
dat <- dat[, -1]

# drop samples with missing response values
dat <- dat[!is.na(dat$response), ]

# drop samples with response values of "Inf" (e.g. GDSC)
mask <- dat$response != Inf
dat <- dat[mask, ]

# drop any features with missing values
dat <- dat[, complete.cases(t(dat))]

# perform feature selection (first attempt)
if (snakemake@config$feature_selection$method == 'boruta') {
  features <- boruta_feature_selection(dat, snakemake) 
} else if (snakemake@config$feature_selection$method == 'rfe') {
  features <- rfe_feature_selection(dat, snakemake) 
} else if (snakemake@config$feature_selection$method == 'none') {
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
    features <- boruta_feature_selection(dat, snakemake) 
  } else if (snakemake@config$feature_selection$fallback == 'rfe') {
    # fallback: rfe
    features <- rfe_feature_selection(dat, snakemake) 
  } else if (snakemake@config$feature_selection$method == 'none') {
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
dat <- dat[, colnames(dat) %in% c('symbol', features, 'response')]

# store result
write_tsv(dat, snakemake@output[[1]])


