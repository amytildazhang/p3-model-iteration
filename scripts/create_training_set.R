#!/bin/env Rscript
#
# combine feature data into a single training set
#
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)

# create a list to store non-empty datasets in
dsets <- list()

data_types <- c('rna', 'cnv', 'variants')

for (i in 1:length(data_types)) {
  # load feature data
  infile <- snakemake@input[[i]]
  dat <- read_tsv(infile, col_type = cols())

  # skip empty datasets
  if (nrow(dat) == 0) {
    print(infile)
    print("nope!")
    next
  }

  # add data type prefix to variables, e.g.
  # "A1BG" -> "cnv_A1BG", "PC1" -> "rna_PC1", etc.
  dat[, 1] <- paste0(data_types[i], '_', as.data.frame(dat)[, 1])
  dsets <- c(dsets, setNames(list(dat), basename(infile)))
}

# get a list of matching sample ids; in most cases (for p3 at least) all feature 
# datasets were created for the same samples;
# one exception to this, however, is the original version of the CNV data from
# the 'max_scores' file, which is missing data for SKMM1_PLB.
shared_sample_ids <- Reduce(intersect, lapply(dsets, function(x) { colnames(x)[-1] }))

print(shared_sample_ids)
print(dsets)

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

# index of response data in snakemake@input
RESPONSE_INPUT_IND <- 4

# add response data and save
response <- read_tsv(snakemake@input[[RESPONSE_INPUT_IND]]) %>%
  filter(cell_line %in% shared_sample_ids)

# response data is stored in an N x 2 dataframe, where the first column lists the
# cell line and the second column the corresponding response values for that cell line
RESPONSE_VALUES_IND <- 2

# get response values for cell lines in the feature data
response <- response[match(shared_sample_ids, response$cell_line), RESPONSE_VALUES_IND]

# add to right side and store result
combined_dat <- cbind(combined_dat, response)
colnames(combined_dat)[ncol(combined_dat)] <- 'response'

write_tsv(combined_dat, snakemake@output[[1]])

