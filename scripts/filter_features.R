#!/bin/env Rscript
#
# Performs basic variance- and correlation-based feature filtering
#
suppressMessages(library(tidyverse))
suppressMessages(library(caret))

options(stringsAsFactors = FALSE)

# load feature data
dat <- read_tsv(snakemake@input[[1]], col_types = cols())


# first column contains the test/train indicator
CV_IND_OFFSET = 1

#
# Separate test/train indicator 
# 
test_index <- CV_IND_OFFSET + which(dat[1, -1] == 0)
cv_mask <- dat[1, ]


# Split into train/test datasets
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
# large feature matrices are split into smaller submatrix blocks 
# for correlation calculations to reduce memory usage and time
if (max_cor < 1 && nrow(dat) > 2) {
  # Use WGCNA::cor as it performed the best in benchmark tests
  suppressMessages(library(WGCNA))

  # get matrix of feature data
  feat_mat <- as.matrix(t(dat[,-1]))

  # large matrices 
  LARGE_MATRIX_WARNING <- 3.5e4

  # max matrix size before it is split into multiple submatrices for correlations
  SPLIT_MATRIX_CUTOFF <- 2.4e4 

  # max size of submatrix blocks
  SUBMATRIX_SIZE <- 1e4

  # maximum number of iterations to split and re-combine matrix into submatrices
  SPLIT_ITER <- 5

  # maximum matrix size to perform correlation calculation on full matrix
  # based on some informal tests correspoding matrix size to memory usage
  MAX_FULL_COR_SIZE <- 3e4

  
  if (ncol(feat_mat) >= LARGE_MATRIX_WARNING) {
    message("Large feature matrix encountered...")
  }

  # initialize vector indicating which columns of feature matrix to keep
  ind <- rep(TRUE, ncol(feat_mat))

  # initialize matrix splitting iterations
  iter <- 0

  # break into chunks until correlation matrix is a manageable size
  while (ncol(feat_mat) > SPLIT_MATRIX_CUTOFF & iter < SPLIT_ITER) {
    print(sprintf("Finding high correlations in submatrices. Matrix size: %d x %d",
                  sum(ind), nrow(feat_mat)))

    # calculate number of submatrices to split full matrix into
    n_waves <- ceiling(ncol(feat_mat) / SUBMATRIX_SIZE)

    # randomly select which columns correspond to which submatrix
    waves <- sample(1:n_waves, size = ncol(feat_mat), replace = TRUE)

    # for each submatrix, calculate correlation matrix and return
    # index of features which have high correlation
    remove <- do.call("c", lapply(1:n_waves, function(i) {
        idx <- waves == i
        cor_mat <- WGCNA::cor(feat_mat[, idx], use = 'pairwise.complete.obs', nThreads = snakemake@threads)

        # return index within feat_mat of features which have high correlation
        which(idx)[caret::findCorrelation(cor_mat, max_cor)]
    }))

    # remove high-correlation features from feat_mat
    feat_mat <- feat_mat[, -remove]

    # mark high-correlation features to discard
    ind[ind][remove] <- FALSE

    iter <- iter + 1
  }

  # if remaining features are of a manageable size, calculate correlaiton matrix for
  # the full feat_mat 
  if (sum(ind) <= MAX_FULL_COR_SIZE) {
    cor_mat <- WGCNA::cor(feat_mat, use = 'pairwise.complete.obs', nThreads = snakemake@threads)
    ind[ind][caret::findCorrelation(cor_mat, max_cor)] <- FALSE
  }

  # remove high-correlation features if any were found
  if (sum(ind) < nrow(dat)) {
    dat <- dat[ind, ]
  }
}


# remove high-correlation features from test dataset and join to train dataset
test_dat <- left_join(dat[, 1], test_dat)

# re-combine with test/train indicator, aligning based on columns
dat <- rbind(cv_mask, cbind(dat, test_dat)[, colnames(cv_mask)])

#
if (max_cor < 1 && nrow(dat) > 2) {
  if (nrow(dat) != sum(ind) + 1) {
    # if duplicated feature names are present, nrow(dat) > sum(ind) + 1
    stop("Lost or added features along the way.")
  } else if (ncol(dat) != length(cv_mask)) {
    # sanity check to ensure the rbind is working as expected
    stop("Lost or added samples along the way.")
  }
}

# store filtered result
write_tsv(dat, snakemake@output[[1]])

