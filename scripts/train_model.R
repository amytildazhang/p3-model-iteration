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
library(caret)
library(parsnip)
library(readr)

set.seed(1)

# load training set
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# DEV
infile <- '../output/gdsc/0.1/08/train/training_sets/selected/GDSC_embelin_ic50_recomputed.tsv.gz'
dat <- read_tsv(infile)

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
}

# store result
write_tsv(mod, snakemake@output[[1]])

sessionInfo()

