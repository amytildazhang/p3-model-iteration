#!/bin/env Rscript
#
# combine feature data into a single training set
#
library(tidyverse)
library(caret)
library(coop)

rna <- read_tsv(snakemake@input[[1]])
cnv <- read_tsv(snakemake@input[[2]])
variants <- read_tsv(snakemake@input[[3]])

# add data type prefix to variables, e.g.
# "A1BG" -> "cnv_A1BG", "PC1" -> "rna_PC1", etc.
rna[, 1] <- paste0('rna_', as.data.frame(rna)[, 1])
cnv[, 1] <- paste0('cnv_', as.data.frame(cnv)[, 1])
variants[, 1] <- paste0('variants_', as.data.frame(variants)[, 1])

# combine into a single dataframe
sample_ids <- intersect(intersect(colnames(rna)[-1], colnames(cnv)[-1]), colnames(variants)[-1])

rna <- rna[, c(colnames(rna)[1], sample_ids)]
cnv <- cnv[, c(colnames(cnv)[1], sample_ids)]
variants <- variants[, c(colnames(variants)[1], sample_ids)]

# remove low variance features from each dataset
row_vars <- apply(rna[, -1], 1, var)
var_cutoff <- quantile(row_vars, snakemake@config$feat_selection_min_var_quantile)
rna <- rna[row_vars >= var_cutoff, ]

row_vars <- apply(cnv[, -1], 1, var)
var_cutoff <- quantile(row_vars, snakemake@config$feat_selection_min_var_quantile)
cnv <- cnv[row_vars >= var_cutoff, ]

row_vars <- apply(variants[, -1], 1, var)
var_cutoff <- quantile(row_vars, snakemake@config$feat_selection_min_var_quantile)
variants <- variants[row_vars >= var_cutoff, ]

# remove correlated features from each dataset;
# perfomring on a per-dataset basis reduces memory requirements significantly and
# should remove most of the correlated features 
cor_mat <- coop::pcor(t(rna[, -1]))
ind <- findCorrelation(cor_mat, snakemake@config$feat_selection_max_cor_rna)
rna <- rna[-ind, ]

cor_mat <- coop::pcor(t(cnv[, -1]))
ind <- findCorrelation(cor_mat, snakemake@config$feat_selection_max_cor_cnv)
cnv <- cnv[-ind, ]

cor_mat <- coop::pcor(t(variants[, -1]))
ind <- findCorrelation(cor_mat, snakemake@config$feat_selection_max_cor_variants)
variants <- variants[-ind, ]

# combine features
dat <- rbind(rna, cnv, variants)

feat_ids <- as.data.frame(dat)[, 1]
dat <- dat[, -1]

dat <- t(dat)
colnames(dat) <- feat_ids

# add response data and save
response <- read_tsv('../data/hmcl/raw/response/NCGC00015038-08_bg_adj_curves_DATA7.tsv.gz') %>%
  filter(cell_line %in% sample_ids)
response <- response[match(sample_ids, response$cell_line), 2]

dat <- cbind(dat, response)
colnames(dat)[ncol(dat)] <- 'response'

write_tsv(dat, snakemake@output[[1]])

