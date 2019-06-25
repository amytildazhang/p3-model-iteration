#!/bin/env Rscript
#
# Perform basic feature selection
#
library(readr)
library(caret)
library(coop)

dat <- read_tsv(snakemake@input[[1]])

# determine input data type from wildcards
if ('rna' %in% names(snakemake@wildcards)) {
  data_type <- 'rna'
} else if ('cnv' %in% names(snakemake@wildcards)) {
  data_type <- 'cnv'
} else if ('var' %in% names(snakemake@wildcards)) {
  data_type <- 'var'
}

# remove low variance features from dataset
row_vars <- apply(dat[, -1], 1, var)
var_cutoff <- quantile(row_vars, snakemake@config$feat_selection[[data_type]][['min_var_quantile']])
dat <- dat[row_vars >= var_cutoff, ]

# remove correlated features from dataset
cor_mat <- coop::pcor(t(dat[, -1]))

ind <- findCorrelation(cor_mat, snakemake@config$feat_selection[[data_type]][['max_cor']])
dat <- dat[-ind, ]

write_tsv(dat, snakemake@output[[1]])

