#!/bin/env Rscript
#
# generate cv fold data subsets
#
suppressMessages(library(tidyverse))
options(stringsAsFactors = FALSE)

# get cv fold indices
cv_fold <- snakemake@params$cv_folds[[snakemake@wildcards$cv]]

# load feature data
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# create test and train folds
train <- dat[cv_fold$train, ]
test <- dat[cv_fold$test, ]

# store result
write_tsv(train, snakemake@output[[1]])
write_tsv(test, snakemake@output[[2]])

