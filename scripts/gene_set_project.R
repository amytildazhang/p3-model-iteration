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

# iterate over gene sets and apply function to data for gene in each set
res <- NULL

for (gset in names(gene_sets)) {
  dat_gset <- dat %>%
    filter(symbol %in% gene_sets[[gset]])

  res <- rbind(res, apply(dat_gset[, -1], 2, snakemake@config$aggregation_func))
}
res <- bind_cols(gene_set = names(gene_sets), as.data.frame(res))

write_tsv(res, snakemake@output[[1]])

