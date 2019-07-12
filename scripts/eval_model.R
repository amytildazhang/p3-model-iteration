#!/bin/env Rscript
#
# Train a model for the specified training set
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
#


source("scripts/model_functions.R")
library(caret)
library(parsnip)
library(tidyverse)

library(rstanarm)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(1)

# load training set
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# identify test/train dataset membership
train_idx <- dat[,2] == 1

# scale and center responses, so that prior scale is less likely to be inappropriate
y <- pull(dat, response)
mean_y <- mean(y[train_idx]); sd_y <- sd(y[train_idx]);
dat$response <- (y - mean_y)/sd_y
y_col <- which(colnames(dat) == "response")

# convert training set feature names to valid column names
cnames <- colnames(dat)
colnames(dat) <- make.names(cnames)

# Choose model to run
if (snakemake@config$model$method == 'random_forest') {
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
     X_test <- cbind(1, as.matrix(dat[!train_idx, -c(1,2,y_col)]))

     if (snakemake@config$model$method == 'linear') {
          if (ncol(dat) > nrow(dat)) stop("Number of features is larger than number of samples--did you mean to run a regularized model?")
          armfit <- stan_glm(response ~ ., dat[train_idx, -c(1:2)], 
                             family = gaussian(),
                             prior = student_t(3, scale = 10),
                             prior_aux = cauchy(0, 10),
                             prior_intercept = student_t(3, scale = 10))

          mod <- draw_post.linear(armfit, dat$sample_id[!train_idx], X_test, pull(dat, y_col)[!train_idx])


     } else if (snakemake@config$model$method == 'bimixture') {
          stdat <- list(n = sum(train_idx),
              y = pull(dat, response)[train_idx],
              X = cbind(1, as.matrix(dat[train_idx, -c(1:2)])),
              p = ncol(dat) -2)
          stanfit <- stan("scripts/bimixture_pointresp.stan",
                data = stdat, chains = 4,
                pars = c("B", "mu", "sigma"),
                control = list(adapt_delta = 0.95,
                               max_treedepth = 20))

          mod <- draw_post.bimodal(stanfit, dat$sample_id[!train_idx], X_test, as.vector(dat[!train_idx, ncol(dat)]))

     }
}

# store result
write_tsv(as_tibble(mod), snakemake@output[[1]])


