#!/bin/env Rscript
#
# Performs basic variance- and correlation-based feature filtering
#
suppressMessages(library(tidyverse))
suppressMessages(library(caret))

options(stringsAsFactors = FALSE)

# load feature data
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

#
# TODO: add comment
#
# do we expect there to be duplicated identifiers at this stage? if so, should we
# perhaps add a step early on to combine or remove such duplicates?.. this way
# we don't have to worry about it throwing off later steps..
#
dat <- dat[!duplicated(pull(dat, 1)), ]

# ID train data to base filtering off of
CV_IND_OFFSET = 1

#
# TODO: add comment
# 
test_index <- CV_IND_OFFSET + which(dat[1, -1] == 0)
cv_indices <- dat[1, ]

test_dat <- dat[-1, c(1, test_index)]

dat <- dat[-1, setdiff(1:ncol(dat), test_index)]

# determine input data type from wildcards ("select_xx_")
data_type <- strsplit(snakemake@rule, '_')[[1]][[2]] 

###########################################
#
# variance-based filtering
#
###########################################
var_quantile <- snakemake@config$feature_filtering[[data_type]][['min_var_quantile']]

# if enabled, remove features with low variance
if (var_quantile > 0 && nrow(dat) > 1) {
  row_vars <- apply(dat[,-1] , 1, var, na.rm = TRUE)
  var_cutoff <- quantile(row_vars, var_quantile)
  dat <- dat[-(which(row_vars < var_cutoff)), ]
}

###########################################
#
# correlation-based filtering
#
###########################################
max_cor <- snakemake@config$feature_filtering[[data_type]][['max_cor']]

# if enabled, remove highly correlated features
if (max_cor < 1 && nrow(dat) > 2) {
  #
  # TODO: explain motivation for using WGCNA::cor
  #
  suppressMessages(library(WGCNA))

  # get matrix of feature data
  feat_mat <- as.matrix(t(dat[,-1]))

  #
  # TODO: comment / rename "magic number" variables
  #
  MAGIC_NUM1 <- 3.5e4
  MAGIC_NUM2 <- 2.4e4
  MAGIC_NUM3 <- 5
  MAGIC_NUM4 <- 3e4
  MAGIC_NUM5 <- 1e4
  
  if (ncol(feat_mat) >= MAGIC_NUM1) {
    message("Large feature matrix encountered...")
  }

  #
  # TODO: comment
  #
  ind <- rep(TRUE, ncol(feat_mat))
  iter <- 0

  # break into chunks until correlation matrix is a manageable size
  while (ncol(feat_mat) > MAGIC_NUM2 & iter < MAGIC_NUM3) {
    print(sprintf("Finding high correlations in submatrices. Matrix size: %d x %d",
                  sum(ind), nrow(feat_mat)))

    #
    # TODO: comment
    #
    n_waves <- ceiling(ncol(feat_mat) / MAGIC_NUM5)
    waves <- sample(1:n_waves, size = ncol(feat_mat), replace = TRUE)

    #
    # TODO: comment
    #
    remove <- do.call("c", lapply(1:n_waves, function(i) {
        idx <- waves == i
        cor_mat <- WGCNA::cor(feat_mat[, idx], use = 'pairwise.complete.obs', nThreads = snakemake@threads)

        #
        # TODO: comment
        #
        which(idx)[caret::findCorrelation(cor_mat, max_cor)]
    }))

    #
    # TODO: comment
    #
    ind[ind][remove] <- FALSE
    feat_mat <- feat_mat[, -remove]

    iter <- iter + 1
  }

  #
  # TODO: comment
  #
  if (sum(ind) <= MAGIC_NUM4) {
    cor_mat <- WGCNA::cor(feat_mat, use = 'pairwise.complete.obs', nThreads = snakemake@threads)
    ind[ind][caret::findCorrelation(cor_mat, max_cor)] <- FALSE
  }

  #
  # TODO: comment
  #
  if (sum(ind) < nrow(dat)) {
    dat <- dat[ind, ]
  }
}

#
# TODO: comment
#
test_dat <- left_join(dat[, 1], test_dat)
dat <- rbind(cv_indices, cbind(dat, test_dat)[, colnames(cv_indices)])

#
# TODO: comment
#
# would we ever expect features/samples to be added at this step?..
#
if (max_cor < 1 && nrow(dat) > 2) {
  if (nrow(dat) != sum(ind) + 1) {
    stop("Lost or added features along the way.")
  } else if (ncol(dat) != ncol(CV_ind)) {
    stop("Lost or added samples along the way.")
  }
}

# store filtered result
write_tsv(dat, snakemake@output[[1]])

