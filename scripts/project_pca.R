#!/bin/env Rscript
#
# Generate a PCA-projected version of feature dataset
#
dat <- read.delim(gzfile(snakemake@input[[1]]), sep = '\t', row.names = 1)

pca <- prcomp(t(dat), scale = snakemake@config$pca_scale)

# the second row of summary(prcomp(..))$importance corresponds to the cumulative
# proportion of variance explained
#
# summary(pca)$importance[1:3, 1:3]
#                             PC1      PC2      PC3
# Standard deviation     5.727381 1.365444 1.105333
# Proportion of Variance 0.762860 0.043360 0.028410
# Cumulative Proportion  0.762860 0.806220 0.834630
VAR_IND <- 3

var_explained <- summary(pca)$importance[VAR_IND, ]

# determine the minimum number of PCs required to explain the specified amount of
# variance in the data
num_pcs <- which(var_explained >= snakemake@config$pca_min_variance)[1]

# get pca-projected version of the data with desired number of PCs
pca_dat <- t(pca$x[, 1:num_pcs])

# save to disk
outfile <- sub('.gz', '', snakemake@output[[1]])
write.table(pca_dat, outfile, sep = '\t', col.names = NA)
system(sprintf("gzip %s", outfile))
