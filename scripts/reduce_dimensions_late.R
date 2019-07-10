#!/bin/env Rscript
#
# Apply dimension reduction to the combined training set
#
options(stringsAsFactors = FALSE)
library(readr)
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

cnames <- colnames(dat)
colnames(dat) <- make.names(cnames)

#
# PCA
#
if (snakemake@config$dimension_reduction_late$method == 'sparse_pls') {
    library(spls)
    coltypes <- sapply(1:ncol(dat), function(col) typeof(dat[[col]]))
    feat <- as.matrix(dat[,coltypes != "character"])
    
    # impute NA and infinite values
    infna_idx <- is.na(feat) | is.infinite(feat)
    infna_cols <- col(feat)[infna_idx]
    infna_rows <- row(feat)[infna_idx]
    if (length(infna_cols) > 0) {
      mean_vals <- sapply(unique(infna_cols), function(j) mean(feat[-infna_rows[infna_cols == j],j], na.rm = T))
      sd_vals <- sapply(unique(infna_cols), function(j) sd(feat[-infna_rows[infna_cols == j],j], na.rm = T))
      na_col_ord <- as.numeric(as.factor(infna_cols))
      feat[infna_idx] <- rnorm(length(na_col_ord), mean = mean_vals[na_col_ord], sd = sd_vals[na_col_ord])
    }

    
    y <- feat[,ncol(feat)]

    # Columns with too many zeros can end up having zero variance in the inner CV, which makes SPLS unhappy
    zero_fvars <- apply(feat[,-ncol(feat)], 2, function(col) sum(col != 0))/nrow(feat) <= 0.12
    if (sum(zero_fvars) > 0)
        message(paste0("Excluding ", sum(zero_fvars), " features with large number of 0s."))

    cv_chs <- cv.spls(feat[, c(!zero_fvars, F)], y,
                               eta = c(seq(0.3, 0.7, by = 0.1), 0.75),  K = c(1:10), fold = 10)
 
    spfit <- spls(feat[, c(!zero_fvars, F)], y, eta = cv_chs$eta.opt, K = cv_chs$K.opt)
    selected <- rownames(spfit$projection)

    latent_feats <- feat[,selected] %*% spfit$projection
    colnames(latent_feats) <- paste0("PC", 1:ncol(latent_feats))
    write_tsv(as.data.frame(cbind(dat[,1], latent_feats)), snakemake@output[[1]])

    mat <- matrix(0, nrow = nrow(feat), ncol = ncol(spfit$projection))
    mat[selected, ] <- spfit$projection
    save(mat, snakemake@output[[2]])
} else {
  stop("Unsupported dimension reduction method specified!")
}
