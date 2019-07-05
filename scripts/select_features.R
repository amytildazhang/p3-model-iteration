#!/bin/env Rscript
#
# Train a random forest model using input training set
#
library(Boruta)
library(readr)

set.seed(1)

# load full training set
dat <- read_tsv(snakemake@input[[1]])

# drop sample ids and convert to a data frame
dat <- as.data.frame(dat)[, -1]

# drop samples with missing response values
dat <- dat[!is.na(dat$response), ]

# drop samples with response values of "Inf" (e.g. GDSC)
mask <- dat$response != Inf
dat <- dat[mask, ]

# drop any features with missing values
dat <- dat[, complete.cases(t(dat))]

RESPONSE_IND <- ncol(dat)

boruta <- Boruta(x = dat[, -RESPONSE_IND], y = dat[, RESPONSE_IND], 
                 num.trees = snakemake@config$feature_selection$num_trees, 
                 mtry = snakemake@config$feature_selection$mtry, doTrace = 2)

selected <- getSelectedAttributes(boruta, withTentative = TRUE)

# remove unselected features
dat <- dat %>%
  select(c('symbol', selected, 'response'))

# store result
write_tsv(dat, snakemake@output[[1]])

