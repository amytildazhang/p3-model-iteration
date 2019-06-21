#!/bin/env Rscript
library(annotables)
library(tidyverse)
options(stringsAsFactors = FALSE)

exclude_cells <- c('KMS21BM_JCRB', 'Karpas417_ECACC', 'OCIMY1_PLB', 'OPM2_DSMZ')

# load raw gene-level (ensembl) counts
dat <- read_tsv('/data/projects/nih/p3/rnaseq/counts_matrix_grch38_featurecounts.tsv')

# drop outlier cell lines
dat <- dat[, !colnames(dat) %in% exclude_cells]

# convert ensembl gene ids -> gene symbols
dat$symbol <- grch38$symbol[match(dat$ensgene, grch38$ensgene)]

#table(is.na(dat$symbol))
# 
# FALSE  TRUE 
# 57861   964 

# remove unmapped genes
dat <- dat[!is.na(dat$symbol), ]

#head(sort(table(dat$symbol), TRUE))
#       Y_RNA Metazoa_SRP          U3      uc_338          U6      snoU13 
#         756         169          51          36          33          31 

# average multi-mapped gene symbols
dat <- dat %>%
  select(-ensgene) %>%
  group_by(symbol) %>%
  summarise_all(funs(mean))

# remove genes with no variance
mask <- apply(dat[, -1], 1, var) > 0

#table(mask)
# mask
# FALSE  TRUE 
#  6028 50206 

dat <- dat[mask, ]

# create a cpm-scaled version of counts
cpm_dat <- dat
cpm_dat[, -1] <- sweep(dat[, -1], 2, colSums(dat[, -1]), '/') * 1E6

write_tsv(dat, '/data/projects/nih/p3-models/data/raw/features/rnaseq/rnaseq_grch38_featurecounts_raw.tsv')
write_tsv(cpm_dat, '/data/projects/nih/p3-models/data/raw/features/rnaseq/rnaseq_grch38_featurecounts_cpm.tsv')

