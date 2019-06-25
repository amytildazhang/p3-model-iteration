#!/bin/env Rscript
#
# aggregate feature data along specified gene sets
#
library(GSEABase)
library(tidyverse)

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

for (gset in names(gene_sets)) {
  dat_gset <- dat %>%
    filter(symbol %in% gene_sets[[gset]])

  # if no genes from gene set were found, continue to next gene set
  if (nrow(dat_gset) == 0) {
    next
  }

  res <- rbind(res, apply(dat_gset[, -1], 2, snakemake@config$aggregation_func))
}
res <- bind_cols(gene_set = names(gene_sets), as.data.frame(res))

write_tsv(res, snakemake@output[[1]])

