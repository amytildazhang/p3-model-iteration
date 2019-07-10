#!/bin/env Rscript
#
# generate cv fold data subsets
#
suppressMessages(library(tidyverse))
options(stringsAsFactors = FALSE)

# get cv fold indices
cv_fold <- snakemake@params$cv_folds[[snakemake@wildcards$drug]][[snakemake@wildcards$cv]]

#
# $train
# [1]  0  1  2  3  4  6 11 13 14 16 19 20 21 22 24 25 26 29 30 32 33 34 36 37 38
# [26] 39 40 41 42
#
# $test
# [1]  5  7  8  9 10 12 15 17 18 23 27 28 31 35
#

# load feature data
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# create test and train folds;

# first column contains identifiers
ID_COL_OFFSET <- 1

# python -> r indexing offset
R_INDEX_OFFSET <- 1

train <- dat[, c(1, cv_fold$train + ID_COL_OFFSET + R_INDEX_OFFSET)]
test <- dat[, c(1, cv_fold$test + ID_COL_OFFSET + R_INDEX_OFFSET)]

# store result
write_tsv(train, snakemake@output[[1]])
write_tsv(test, snakemake@output[[2]])

