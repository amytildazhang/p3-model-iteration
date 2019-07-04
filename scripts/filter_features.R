#!/bin/env Rscript
#
# Perform basic feature selection
#
suppressMessages(library(readr))
suppressMessages(library(caret))
suppressMessages(library(coop))

options(stringsAsFactors = FALSE)

dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# determine input data type from wildcards ("select_xx_")
data_type <- strsplit(snakemake@rule, '_')[[1]][[2]] 

# remove low variance features from dataset
row_vars <- apply(dat[, -1], 1, var, na.rm = TRUE)

var_cutoff <- quantile(row_vars, snakemake@config$feat_selection[[data_type]][['min_var_quantile']])
dat <- dat[row_vars >= var_cutoff, ]

# remove correlated features from dataset
cor_mat <- coop::pcor(t(dat[, -1]), use = 'pairwise.complete')

ind <- findCorrelation(cor_mat, snakemake@config$feat_selection[[data_type]][['max_cor']])

if (length(ind) > 0) {
  dat <- dat[-ind, ]
}

write_tsv(dat, snakemake@output[[1]])
