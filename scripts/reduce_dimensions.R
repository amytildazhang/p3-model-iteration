#!/bin/env Rscript
#
# Apply dimension reduction to the combined training set
#
suppressMessages(library(tidyverse))
options(stringsAsFactors = FALSE)

dat <- read_tsv(snakemake@input[[1]], col_types = cols())


# separate out cross-validation indicator and identify the train set
train_idx <- pull(dat, train_idx) == 1
y <- pull(dat, response)
y_col <- which(colnames(dat) == "response")
feat <- as.matrix(dat[,-c(1,2,y_col)])


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

METHOD <- snakemake@wildcards$dim_red

#
# Sparse PLS
#
if (METHOD == 'sparse_pls') {
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
    colnames(latent_feats) <- paste0("LF", 1:ncol(latent_feats))

    # save projection matrix
    projection <- t(spfit$projection)
    colnames(projection) <- selected 
    write_tsv(as.data.frame(projection), snakemake@output[[2]])
} else if (METHOD == 'mofa') {
    library(MOFA)

    # MoFA prefers having sample names
    rownames(feat) <- pull(dat, sample_id) 

    # separate data into view-specific matrices
    view_labs <- str_split(colnames(feat), "_", simplify = T)[,1]   
    views <- unique(view_labs)    
    view_mats <- map(views, function(vw) {
        mat <- t(feat[, view_labs == vw]) # must fit to full data
        rownames(mat) <- NULL
        mat
    }) %>% set_names(views)

    # create MOFA object
    MOFAobj <- createMOFAobject(view_mats)

    # pull MoFA settings from the config file
    mofa_sets <- snakemake@config$dimension_reduction_late$mofa
    
    # MoFA settings: number of factors, thresholds, iterations, seed
    ModelOptions <- getDefaultModelOptions(MOFAobj)
    for (opt in names(mofa_sets$modeloptions)) ModelOptions[[opt]] <- mofa_sets$modeloptions[[opt]]
    ModelOptions$numFactors <- min(ceiling(nrow(feat)/2), ModelOptions$numFactors)    
 
    TrainOptions <- getDefaultTrainOptions()
    for (opt in names(mofa_sets$trainoptions)) TrainOptions[[opt]] <- mofa_sets$trainoptions[[opt]]

    # run MoFA
    MOFAobj <- prepareMOFA(MOFAobj, ModelOptions = ModelOptions, TrainOptions = TrainOptions)
    MOFAobj <- runMOFA(MOFAobj)

    # use latent factors as features
    latent_feats <- getFactors(MOFAobj)

    # save weight matrices (factor x feature)
    weights <- getWeights(MOFAobj)
    weights <- lapply(views, function(vw) {
        mat <- weights[[vw]]
        rownames(mat) <- colnames(feat)[view_labs == vw]
        t(mat)
    })
    write_tsv(as.data.frame(do.call("cbind", weights)), snakemake@output[[2]])

} else {
  stop("Unsupported dimension reduction method specified!")
}

# save latent features (sample x factor)
write_tsv(as.data.frame(cbind(dat[,c(1,2,ncol(dat))], latent_feats)), snakemake@output[[1]])
