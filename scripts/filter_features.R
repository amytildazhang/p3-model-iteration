#!/bin/env Rscript
#
# Perform variance- and correlated-based feature filtering
#
suppressMessages(library(tidyverse))
suppressMessages(library(caret))
suppressMessages(library(coop))

options(stringsAsFactors = FALSE)

dat <- read_tsv(snakemake@input[[1]], col_types = cols())


# ID train data to base filtering off of
CV_IND_OFFSET = 1
test_idx <- CV_IND_OFFSET + which( dat[1, -1] == 0)
CV_ind <- dat[1,]
test_dat <- dat[-1, c(1, test_idx)]
dat <- dat[-1, setdiff(1:ncol(dat), test_idx)]


# determine input data type from wildcards ("select_xx_")
data_type <- strsplit(snakemake@rule, '_')[[1]][[2]] 

# remove low variance features from dataset
var_quantile <- snakemake@config$feature_filtering[[data_type]][['min_var_quantile']]

if (var_quantile > 0 && nrow(dat) > 1) {
  row_vars <- apply(dat[,-1] , 1, var, na.rm = TRUE)
  var_cutoff <- quantile(row_vars, var_quantile)
  dat <- dat[-(which(row_vars < var_cutoff)), ]
}

# remove correlated features from dataset
max_cor <- snakemake@config$feature_filtering[[data_type]][['max_cor']]

if (max_cor < 1 && nrow(dat) > 2) {
  cor_mat <- coop::pcor(t(dat[, -1]), use = 'pairwise.complete')

  ind <- findCorrelation(cor_mat, max_cor)
  if (length(ind) > 0) {
    dat <- dat[-ind, ]
  }
}

dat <- bind_rows(CV_ind, left_join(dat, test_dat, by = "gene_set"))

write_tsv(dat, snakemake@output[[1]])

