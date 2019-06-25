#!/bin/env Rscript
#
# aggregate feature data along specified gene sets
#
library(GSEABase)
library(tidyverse)

# determine input data type from wildcards
if ('rna' %in% names(snakemake@wildcards)) {
  data_type <- 'rna'
} else if ('cnv' %in% names(snakemake@wildcards)) {
  data_type <- 'cnv'
} else if ('var' %in% names(snakemake@wildcards)) {
  data_type <- 'var'
}

# load feature data
dat <- read_tsv(snakemake@input$features)

gene_sets <- geneIds(getGmt(gzfile(snakemake@input$gene_set)))

# remove gene set :length suffixes, if present
names(gene_sets) <- sub(':\\d+$', '', names(gene_sets)) 

# remove gene set weights, if present
# e.g. "ANXA1,1.0" -> "ANXA1"
gene_sets <- lapply(gene_sets, function(x) { sub(',\\d+\\.\\d+$', '', x) }) 

# exclude any gene sets with fewer than the required number of genes
set_sizes <- lapply(gene_sets, length)
mask <- set_sizes >= snakemake@config$gene_set_min_size

gene_sets <- gene_sets[mask]

# iterate over gene sets and apply function to data for gene in each set
res <- NULL
gset_names <- c()

# determine which aggregation function to use
agg_func <- snakemake@config$aggregation_funcs[[data_type]]

for (gset in names(gene_sets)) {
  dat_gset <- dat %>%
    filter(symbol %in% gene_sets[[gset]])

  # if no genes from gene set were found, continue to next gene set
  if (nrow(dat_gset) == 0) {
    next
  }
  gset_names <- c(gset_names, gset)

  res <- rbind(res, apply(dat_gset[, -1], 2, agg_func))
}

# drop any rows with zero variance (uninformative)
row_vars <- apply(res, 1, var)
res <- res[row_vars != 0, ]

res <- bind_cols(gene_set = gset_names, as.data.frame(res))

write_tsv(res, snakemake@output[[1]])

