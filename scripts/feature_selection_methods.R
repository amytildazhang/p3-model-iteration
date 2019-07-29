#
# Feature selection methods
#
# - Boruta
# - RFE/RF
# - Distance correlation

#
# Boruta
#
boruta_feature_selection <- function(dat, snakemake) {
  library(Boruta)

  RESPONSE_IND <- ncol(dat)

  # perform boruta feature selection
  boruta <- Boruta(x = dat[, -RESPONSE_IND], y = dat[, RESPONSE_IND], 
                  num.trees = snakemake@config$feature_selection$num_trees, 
                  mtry = snakemake@config$feature_selection$mtry, doTrace = 1,
                  num.threads = snakemake@config$num_threads$feature_selection)

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
  library(dplyr)
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
  MAX_THREADS <- max(1, min(detectCores() - 4, snakemake@config$num_threads$feature_selection, na.rm = TRUE))

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

  # if too few variables found, use ranking to retrieve top N vars intead
  if (num_features < snakemake@config$feature_selection$min_features) {
    features <- rfe_res$variables %>%
      group_by(var) %>%
      summarize(mean_score = mean(Overall)) %>%
      arrange(desc(mean_score)) %>%
      head(snakemake@config$feature_selection$min_features) %>%
      pull(var)
  } else {
    # otherwise, return selected features
    features <- head(rfe_res$optVariables, num_features)
  }

  features
}


#
# Distance correlation
#

#'
#' Returns the p/log(p) features most correlated with the response, based on
#' distance correlation. 
#' --------
#'
#' Distance correlation (DC) is 0 only if two random variables are independent and is based on 
#' Euclidean distances between samples rather than moments (Szekely, Rizzo, Bakirov 2007). It is
#' also a strictly increasing function of Pearson correlation.
#' 
#' Li, Zhong, Zhu (2012) prove all features in the true model have DC greater than some threshold 
#' \code{c}, but \code{c} is an unknown value based on sample and feature size. A suggested alternative is to choose the top \code{n}
#' features based on distance correlation.Their DC sure independence screening (DC-SIS) method
#' can handle both multivariate response and grouped predictors and is a model-free SIS procedure.
#'
#' 
#' Parameters
#' ----------
#'
#' @param dat A dataframe containing features and response in the last column.
#' @param snakemake \code{snakemake} object containing number of threads for parallelization.
#'
#' Returns
#' -------
#'
#' @return A character vector of feature names with the highest distance correlations to the response. 
#'
#' References
#' ----------
#' 1. Székely, Gábor J., Maria L. Rizzo, and Nail K. Bakirov. "Measuring and testing dependence
#' by correlation of distances." The Annals of Statistics 35.6 (2007): 2769-2794.
#' 2. Li, Runze, Wei Zhong, and Liping Zhu. "Feature screening via distance correlation learning."
#' Journal of the American Statistical Association 107.499 (2012): 1129-1139. #'
dcor_feature_selection <- function(dat, snakemake) {
  library(snow) # manages parallelization

  # pull out response
  Y <- dat[, ncol(dat)]

  # create parallelization cluster instance
  cl <- makeCluster(snakemake@threads)

  # load library and export objects to worker nodes
  clusterCall(cl, function() { library(energy) })
  clusterExport(cl, list("Y", "dat"))

  # start and end indices within dat to calculate distance correlation for
  FEAT_START_IND <- 1
  FEAT_END_IND <- ncol(dat) - 1

  # numeric vector of values [0, 1], where 0 indicates independence between
  # two random variables
  dist_cors <- parSapply(cl, FEAT_START_IND:FEAT_END_IND, function(j) {
    dcor(pull(dat, j), Y)
  })

  # 
  stopCluster(cl)

  # set lower bound for DC based on the DC for the p/log(p)th feature
  num_feats <- length(dist_cors)
  min_val <- abs(sort(-dist_cors)[round(num_feats / log(num_feats))])

  features <- colnames(dat)[dist_cors >= min_val]

  features
}