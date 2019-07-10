#!/bin/env Rscript
#
# Apply dimension reduction to each individual dataset
#
suppressMessages(library(tidyverse))
options(stringsAsFactors = FALSE)

dat <- read_tsv(snakemake@input[[1]], col_types = cols())

#
# PCA
#
if (snakemake@config$dimension_reduction_early$method == 'pca') {
  pca <- prcomp(t(dat[, -1]), scale = snakemake@config$dimension_reduction_early$scale)

  # the second row of summary(prcomp(..))$importance corresponds to the cumulative
  # proportion of variance explained
  #
  # summary(pca)$importance[1:3, 1:3]
  #                             PC1      PC2      PC3
  # Standard deviation     5.727381 1.365444 1.105333
  # Proportion of Variance 0.762860 0.043360 0.028410
  # Cumulative Proportion  0.762860 0.806220 0.834630
  VAR_IND <- 3

  var_explained <- summary(pca)$importance[VAR_IND, ]

  # determine the minimum number of PCs required to explain the specified amount of
  # variance in the data
  num_pcs <- which(var_explained >= snakemake@config$dimension_reduction_early$min_variance)[1]

  # get pca-projected version of the data with desired number of PCs
  pca_dat <- t(pca$x[, 1:num_pcs])

  # save to disk
  pca_dat %>%
    as.data.frame %>%
    rownames_to_column('PC') %>%
    write_tsv(snakemake@output[[1]])
} else {
  stop("Unsupported dimension reduction method specified!")
}
