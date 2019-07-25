#!/bin/env Rscript
#
# Perform variance- and correlated-based feature filtering
#
suppressMessages(library(tidyverse))
suppressMessages(library(caret))

options(stringsAsFactors = FALSE)

dat <- read_tsv(snakemake@input[[1]], col_types = cols())
dat <- dat[!duplicated(pull(dat, 1)), ]

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
    fmat <- as.matrix(t(dat[,-1]))
    if (ncol(fmat) >= 3.5e4) message("Large feature matrix")
  
    ind <- rep(T, ncol(fmat))
    iter <- 0

    # break into chunks until correlation matrix is a manageable size
    while (ncol(fmat) > 2.4e4 & iter < 5) {
        print(paste("Finding high correlations in submatrices. Matrix size", sum(ind), "x", nrow(fmat)))
        n_waves <- ceiling(ncol(fmat)/1e4)
        waves <- sample(1:n_waves, size = ncol(fmat), replace = T)

        remove <- do.call("c", lapply(1:n_waves, function(i) {
            idx <- waves == i
            cor_mat <- WGCNA::cor(fmat[, idx], use = 'pairwise.complete.obs', nThreads = snakemake@threads)
            which(idx)[caret::findCorrelation(cor_mat, max_cor)]
        }))
        ind[ind][remove] <- F
        fmat <- fmat[, -remove]
        iter <- iter + 1

    }

  
    if (sum(ind) <= 3e4) {
        cor_mat <- WGCNA::cor(fmat, use = 'pairwise.complete.obs', nThreads = snakemake@threads)
        ind[ind][caret::findCorrelation(cor_mat, max_cor)] <- F
    }
    if (sum(ind) < nrow(dat)) {
        dat <- dat[ind, ]
    }
}

test_dat <- left_join(dat[,1], test_dat)
dat <- rbind(CV_ind, cbind(dat, test_dat)[, colnames(CV_ind)])


if (max_cor < 1 && nrow(dat) > 2) {
    if (nrow(dat) != sum(ind) + 1) stop("Lost or added features along the way.")
    if (ncol(dat) != ncol(CV_ind)) stop("Lost or added samples along the way.")
}

write_tsv(dat, snakemake@output[[1]])

