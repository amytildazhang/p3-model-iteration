#!/bin/env Rscript
#
# Apply dimension reduction to the combined training set
#
options(stringsAsFactors = FALSE)

dat <- read_tsv(snakemake@input[[1]], col_types = cols())

#
# PCA
#
if (snakemake@config$dimension_reduction_late$method == 'pls') {
  #
  # TODO...
  #
} else {
  stop("Unsupported dimension reduction method specified!")
}
