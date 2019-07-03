#!/bin/env Rscript
#
# combine feature data into a single training set
#
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)

# create a list to store non-empty datasets in
dsets <- list()

# add data-type prefixes (e.g. "rna_") to feature names
data_types <- c('rna', 'cnv', 'var')

for (data_type in data_types) {
  # load feature data
  infile <- snakemake@input[[data_type]]
  dat <- read_tsv(infile, col_type = cols())

  # skip empty datasets
  if (nrow(dat) == 0) {
    print(sprintf("Skipping %s: no features remaining after filtering...", data_type))
    next
  }

  # add data type prefix to variables, e.g.
  # "A1BG" -> "cnv_A1BG", "PC1" -> "rna_PC1", etc.
  dat[, 1] <- paste0(data_type, '_', as.data.frame(dat)[, 1])
  dsets <- c(dsets, setNames(list(dat), basename(infile)))
}

# get a list of matching sample ids; in most cases (for p3 at least) all feature 
# datasets were created for the same samples;
# one exception to this, however, is the original version of the CNV data from
# the 'max_scores' file, which is missing data for SKMM1_PLB.
shared_sample_ids <- Reduce(intersect, lapply(dsets, function(x) { colnames(x)[-1] }))

# normalize columns across feature datasets
for (dset in names(dsets)) {
  id_col <- colnames(dsets[[dset]])[1]
  dsets[[dset]] <- dsets[[dset]][, c(id_col, shared_sample_ids)]
}

# combine features
combined_dat <- do.call(rbind, dsets)

# transpose data and set column names to feature ids
feat_ids <- as.data.frame(combined_dat)[, 1]
combined_dat <- combined_dat[, -1]

combined_dat <- t(combined_dat)
colnames(combined_dat) <- feat_ids

# add response data and save
response <- read_tsv(snakemake@input[['response']], col_type = cols()) %>%
  select(shared_sample_ids)

response <- as.vector(t(as.data.frame(response)[1, ]))

# add response and sample id columns and store result
combined_dat <- as.data.frame(cbind(sample_id = shared_sample_ids, cbind(combined_dat, response)))
write_tsv(combined_dat, snakemake@output[[1]])

