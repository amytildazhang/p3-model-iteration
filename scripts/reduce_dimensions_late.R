#!/bin/env Rscript
#
# Apply dimension reduction to the combined training set
#
suppressMessages(library(tidyverse))
options(stringsAsFactors = FALSE)

dat <- read_tsv(snakemake@input[[1]], col_types = cols())


# separate out cross-validation indicator and identify the train set
train_idx <- pull(dat, 2) == 1
feat <- as.matrix(dat[,3:(ncol(dat)-1)])
y <- pull(dat, ncol(dat))


# impute NA and infinite values
# Unfortunately there is no "infinite.rm = T" option, so ID those manually
infna_idx <- is.na(feat) | is.infinite(feat)
infna_cols <- col(feat)[infna_idx]
infna_rows <- row(feat)[infna_idx]

if (length(infna_cols) > 0) {
  # get mean and sd of non NA/infinite values in each column
  mean_vals <- sapply(unique(infna_cols), function(j) mean(feat[-infna_rows[infna_cols == j],j], na.rm = T))
  sd_vals <- sapply(unique(infna_cols), function(j) sd(feat[-infna_rows[infna_cols == j],j], na.rm = T))

  na_col_ord <- as.numeric(as.factor(infna_cols))
  feat[infna_idx] <- rnorm(length(na_col_ord), mean = mean_vals[na_col_ord], sd = sd_vals[na_col_ord])
}


#
# Sparse PLS
#
if (snakemake@config$dimension_reduction_late$method == 'sparse_pls') {
    library(spls)
    set.seed(2897)
    # columns with too many zeros can end up having zero variance in the inner CV, which makes SPLS unhappy
    zero_fvars <- apply(feat[train_idx,], 2, function(col) sum(col != 0))/sum(train_idx) <= 0.15
    if (sum(zero_fvars) > 0) {
        message(paste0("Excluding ", sum(zero_fvars), " features with large number of 0s."))
	feat <- feat[, !zero_fvars]
    }
    # Select number of components and sparsity penalization (eta) by cross-validation
    cv_chs <- cv.spls(feat[train_idx, ], y[train_idx],
                               eta = c(seq(0.3, 0.7, by = 0.1), 0.75),  K = c(1:10), fold = 10)

    # Fit sparse PLS using selected number of components and eta
    spfit <- spls(feat[train_idx, ], y[train_idx], eta = cv_chs$eta.opt, K = cv_chs$K.opt)
    
    # Project entire feature matrix onto principal components
    selected <- rownames(spfit$projection)
    latent_feats <- feat[,selected] %*% spfit$projection
    colnames(latent_feats) <- paste0("PC", 1:ncol(latent_feats))

    # save
    write_tsv(as.data.frame(cbind(dat[,c(1,2,ncol(dat))], latent_feats)), snakemake@output[[1]])
} else {
  stop("Unsupported dimension reduction method specified!")
}
