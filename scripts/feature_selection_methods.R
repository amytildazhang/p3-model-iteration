#
# Feature selection methods
#
# - Boruta
# - RFE/RF
#

#
# Boruta
#
boruta_feature_selection <- function(dat, snakemake) {
  library(Boruta)

  RESPONSE_IND <- ncol(dat)

  # perform boruta feature selection
  boruta <- Boruta(x = dat[, -RESPONSE_IND], y = dat[, RESPONSE_IND], 
                  num.trees = snakemake@config$feature_selection$num_trees, 
                  mtry = snakemake@config$feature_selection$mtry, doTrace = 1)

  # if development mode is enabled, save model
  if (snakemake@config$debug) {
    save(boruta, file = sub('.tsv.gz', '.rda', snakemake@output[[1]]))
  }

  # get selected features
  features <- getSelectedAttributes(boruta, withTentative = TRUE)

  num_features <- min(length(features), snakemake@config$feature_selection$max_features)

  head(features, num_features)
}

#
# RFE
#
rfe_feature_selection <- function(dat, snakemake) {

  library(caret)
  library(doParallel)
  library(randomForest)

  # use the same cross validation settings as the main outer CV loop
  num_folds   <- snakemake@config$cross_validation$num_splits
  num_repeats <- snakemake@config$cross_validation$num_repeats

  # number of features to include at each iteration
  #num_features <- round(c(ncol(dat) / seq(2, 10, by = 2),
  #                        rev(seq(1, ncol(dat) / 10, by = 100))))

  num_features <- seq(1, round((ncol(dat) - 1)^(1/4)))^(4)

  # for a dataset with 10,000 features, this results in an rfe sizes vector of
  # [1]     1    16    81   256   625  1296  2401  4096  6561 10000


  # generate random seeds
  seeds <- vector(mode = "list", length = num_folds + 1)

  for(i in 1:num_folds) {
    seeds[[i]] <- sample.int(1E6, length(num_features) + 1)
  }
  seeds[[num_folds + 1]] <- sample.int(1E6, 1)

  # rfe settings
  control <- rfeControl(functions = rfFuncs, method = "cv", number = num_folds, 
                        repeats = num_repeats, seeds = seeds, verbose = TRUE)

  # register parallel back-end
  MAX_THREADS <- 16
  MAX_THREADS <- max(1, min(detectCores() - 4, MAX_THREADS, na.rm = TRUE))

  # create paralellization handler
  cl <- makeCluster(MAX_THREADS, outfile='')
  registerDoParallel(cl)

  RESPONSE_IND <- ncol(dat)

  rfe_res <- rfe(dat[, -RESPONSE_IND], dat[, RESPONSE_IND], sizes = num_features, 
                 rfeControl = control)
  stopCluster(cl)

  # if development mode is enabled, save model
  if (snakemake@config$debug) {
    save(rfe_res, file = sub('.tsv.gz', '.rda', snakemake@output[[1]]))
  }

  # determine number of features to keep
  num_features <- min(length(rfe_res$optVariables), snakemake@config$feature_selection$max_features)

  head(rfe_res$optVariables, num_features)
}
