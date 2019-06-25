#!/bin/env Rscript
#
# Perform basic feature selection
#
library(readr)
library(caret)
library(coop)

dat <- read_tsv(snakemake@input[[1]])

# remove low variance features from dataset
row_vars <- apply(dat[, -1], 1, var)
var_cutoff <- quantile(row_vars, snakemake@config$feat_selection_min_var_quantile)
dat <- dat[row_vars >= var_cutoff, ]

# remove correlated features from dataset
cor_mat <- coop::pcor(t(dat[, -1]))

# feature-type specific filtering (not implemented; can add later as-needed)
#feat_type <- snakemake@params$feature_type
#ind <- findCorrelation(cor_mat, snakemake@config$feat_selection_max_cor[[feat_type]])

ind <- findCorrelation(cor_mat, snakemake@config$feat_selection_max_cor)
dat <- dat[-ind, ]

write_tsv(dat, snakemake@output[[1]])

