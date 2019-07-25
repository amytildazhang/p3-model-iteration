#!/bin/env Rscript
#
# generate cv fold data subsets
#
suppressMessages(library(tidyverse))
options(stringsAsFactors = FALSE)

# load feature data
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# get cv fold indices
if (snakemake@wildcards$cv == "alldata") {
  cv_fold <- list(train = 0:ncol(dat))
} else {
  cv_fold <- snakemake@params$cv_folds[[snakemake@wildcards$drug]][[snakemake@wildcards$cv]]
}
#
# $train
# [1]  0  1  2  3  4  6 11 13 14 16 19 20 21 22 24 25 26 29 30 32 33 34 36 37 38
# [26] 39 40 41 42
#
# $test
# [1]  5  7  8  9 10 12 15 17 18 23 27 28 31 35
#

# first column contains identifiers
ID_COL_OFFSET <- 1

# python -> r indexing offset
R_INDEX_OFFSET <- 1


# create indicator variable for test and train folds;
train_ind <- as.numeric(2:ncol(dat) %in% (cv_fold$train + ID_COL_OFFSET + R_INDEX_OFFSET))

# change ID column to include name for indicator var
feat_names <- c("train_idx", pull(dat,1 ))

# save indicator var as first row of data frame
dat_cv <- cbind(symbol = feat_names, as.data.frame(rbind(train_ind, as.matrix(dat[,-1]))))
colnames(dat_cv) <- colnames(dat)


# store result
write_tsv(dat_cv, snakemake@output[[1]])

