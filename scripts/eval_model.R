#!/bin/env Rscript
################################################################################
#
# Trains a model for the specified training set
#
# In order to provide consistent output from various models, the tidymodels package is
# used.
#
# Models currently supported by tidymodels:
#
#  https://tidymodels.github.io/parsnip/articles/articles/Models.html 
#
# Information on addinng support for custom methods:
#
#  https://tidymodels.github.io/parsnip/articles/articles/Scratch.html
#  https://tidymodels.github.io/model-implementation-principles/
#
################################################################################

# if running on Biowulf, load openblas + gcc/7.3.0
if (startsWith(Sys.info()[['nodename']], 'cn')) {
  system("module load openblas gcc/7.3.0")
}

source("scripts/model_functions.R")

library(caret)
library(parsnip)
library(tidyverse)
library(rstanarm)

set.seed(1)

# load training set
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# column containing sample names
SAMPLE_IND <- 1


# identify test/train dataset membership
CV_IND <- 2 # column with test/train indicator
train_idx <- dat[, CV_IND] == 1

# identify whether we're fitting to all data or not
FIT_FULL <- all(train_idx)

# scale and center responses, so that prior scale is less likely to be inappropriate
y <- pull(dat, response)
mean_y <- mean(y[train_idx]); sd_y <- sd(y[train_idx]);
dat$response <- (y - mean_y)/sd_y
y_col <- which(colnames(dat) == "response")

# convert training set feature names to valid column names
cnames <- colnames(dat)
colnames(dat) <- make.names(cnames)

MODEL <- snakemake@wildcards$model

# Choose model to run
if (MODEL == 'random_forest') {
  # use caret to tune model hyperparameters
	train_control <- trainControl(method = "cv", number = 10, savePred = TRUE,
															  classProb = FALSE)

  caret_mod <- train(response ~ ., data = dat, method = "ranger", metric = 'Rsquared',
                     trControl = train_control, importance = 'permutation',
										 num.trees = snakemake@config$model$num_trees,
                     num.threads = snakemake@config$num_threads$train_model)

  # next, generate a tidy version of the tuned model using parsnip
	mod <- rand_forest(mode = "regression", mtry = caret_mod$bestTune$mtry, 
										 trees = snakemake@config$model$num_trees, 
									   min_n = caret_mod$bestTune$min.node.size) %>%
		set_engine("ranger", importance = 'permutation', seed = 1,
							 splitrule = caret$bestTune$splitrule,
               num.threads = snakemake@config$num_threads$train_model) %>%
		fit(response ~ ., data = dat)
} else {

  
  if (!FIT_FULL) {
    # create design matrix for test data; use to calculate model evaluation metrics
    X_test <- cbind(1, as.matrix(dat[!train_idx, -c(SAMPLE_IND, CV_IND, y_col)]))
  }

  # dimension of features (includes intercept)
  p <- ncol(dat) - 2

  if (MODEL == 'linear') {
    if (ncol(dat) > nrow(dat)) {
      stop("Number of features is larger than number of samples--did you mean to run a regularized model?")
    }

    # fit simple linear model using `rstanarm` package
    armfit <- stan_glm(response ~ ., dat[train_idx, -c(SAMPLE_IND, CV_IND)], 
                      # define prior densities
                      family = gaussian(),
                      prior = student_t(3, scale = 10),
                      prior_aux = cauchy(0, 10),
                      prior_intercept = student_t(3, scale = 10))

    if (FIT_FULL) {
      # save stanfit object
      mod <- armfit$stanfit
    } else {
      # save posterior predictive samples and log(CPO)
      mod <- draw_post.linear(armfit, dat$sample_id[!train_idx], X_test, pull(dat, y_col)[!train_idx])
    }
  } else if (MODEL == 'bimixture') {
    library(rstan)
    options(mc.cores = parallel::detectCores())

    rstan_options(auto_write = TRUE)

    # create data object to pass to Stan
    stdat <- list(n = sum(train_idx),
                  y = pull(dat, response)[train_idx],
                  X = cbind(1, as.matrix(dat[train_idx, -c(SAMPLE_IND, CV_IND, y_col)])),
                  p = p)        

    # fit model using defined script; prior/posterior densities defined in script
    stanfit <- stan("scripts/bimixture_pointresp.stan",
                    data = stdat, chains = 4,
                    pars = c("B", "mu", "sigma"), 
                    control = list(adapt_delta = 0.95, max_treedepth = 20))

    if (FIT_FULL) {
      # save stanfit object
      mod <- stanfit
    } else {
      # save posterior predictive samples and log(CPO)
      mod  <- draw_post.bimodal(stanfit, dat$sample_id[!train_idx], X_test, deframe(dat[!train_idx, y_col]))
    }
  } else if (MODEL == 'horseshoe') {
    # recommend for unsupervised dimension reduction methods 
    # places a Finnish horseshoe prior on all features
    
    EFF_FEATURES_GUESS <- 4  # estimate for number of relevant covariates
    p0 <- ifelse(p < EFF_FEATURES_GUESS, ceiling(p/2), EFF_FEATURES_GUESS)  
    # prior parameter based on estimate for number of relevant covarites
    tau0 <- p0/(p - p0) * 1/sqrt(sum(train_idx)) # 
    

    # fit model
    armfit <- stan_glm(response ~ ., dat[train_idx, -c(SAMPLE_IND, CV_IND)],
                      # define prior densities
                       family = gaussian(),
                       prior = hs(df=1, global_df=1, global_scale = tau0),
                       prior_aux = cauchy(0, 10),
                       prior_intercept = student_t(3, scale = 10))

    if (FIT_FULL) {
      # save stanfit object
      mod <- armfit$stanfit
    } else {
      # save posterior predictive samples and log(CPO)
      mod <- draw_post.linear(armfit, dat$sample_id[!train_idx], X_test, pull(dat, y_col)[!train_idx])
    }
  }
}

# store result, either as stanfit object or posterior predictive samples
if (FIT_FULL) {
  saveRDS(mod, snakemake@output[[1]])
} else {
  write_tsv(as_tibble(mod), snakemake@output[[1]]) 
}

