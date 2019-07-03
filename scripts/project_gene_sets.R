#!/bin/env Rscript
#
# aggregate feature data along specified gene sets
#
suppressMessages(library(GSEABase))
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)

# determine input data type from wildcards ("create_rna_gene_sets")
data_type <- strsplit(snakemake@rule, '_')[[1]][[2]] 

# load feature data
dat <- read_tsv(snakemake@input[[1]], col_types = cols())

# load gene sets
gene_sets = c()

message(sprintf("Loading gene sets for %s", data_type))

for (infile in Sys.glob('gene_sets/*.gmt.gz')) {
  fp <- gzfile(infile)
  gene_sets <- c(gene_sets, geneIds(getGmt(fp)))
  close(fp)
}

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

message(sprintf('Aggregating %s along gene sets...', data_type))

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
mask <- apply(res, 1, var) > 0
res <- res[mask, ]
gset_names <- gset_names[mask]

message(sprintf('Saving gene sets aggregated %s data...', data_type))

res <- bind_cols(gene_set = gset_names, as.data.frame(res))

write_tsv(res, snakemake@output[[1]])

