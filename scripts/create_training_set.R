#!/bin/env Rscript
#
# combine feature data into a single training set
#
library(tidyverse)
library(caret)
library(coop)

rna <- read_tsv(snakemake@input[[1]])
cnv <- read_tsv(snakemake@input[[2]])
variants <- read_tsv(snakemake@input[[3]])

# add data type prefix to variables, e.g.
# "A1BG" -> "cnv_A1BG", "PC1" -> "rna_PC1", etc.
rna[, 1] <- paste0('rna_', as.data.frame(rna)[, 1])
cnv[, 1] <- paste0('cnv_', as.data.frame(cnv)[, 1])
variants[, 1] <- paste0('variants_', as.data.frame(variants)[, 1])

# combine into a single dataframe
sample_ids <- intersect(intersect(colnames(rna)[-1], colnames(cnv)[-1]), colnames(variants)[-1])

rna <- rna[, c(colnames(rna)[1], sample_ids)]
cnv <- cnv[, c(colnames(cnv)[1], sample_ids)]
variants <- variants[, c(colnames(variants)[1], sample_ids)]

# combine features
dat <- rbind(rna, cnv, variants)

feat_ids <- as.data.frame(dat)[, 1]
dat <- dat[, -1]

dat <- t(dat)
colnames(dat) <- feat_ids

# add response data and save
response <- read_tsv(snakemake@input[[4]]) %>%
  filter(cell_line %in% sample_ids)

RESPONSE_VALUES_IND <- 2
response <- response[match(sample_ids, response$cell_line), RESPONSE_VALUES_IND]

dat <- cbind(dat, response)
colnames(dat)[ncol(dat)] <- 'response'

write_tsv(dat, snakemake@output[[1]])

