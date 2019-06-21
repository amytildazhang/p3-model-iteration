#!/bin/env Rscript
library(annotables)
library(tidyverse)
options(stringsAsFactors = FALSE)

exclude_cells <- c('KMS21BM_JCRB', 'Karpas417_ECACC', 'OCIMY1_PLB', 'OPM2_DSMZ')

# load raw gene-level (ensembl) counts
dat <- read_tsv('/data/projects/nih/p3/cnv/cnv_gene_max_scores.tab')

# drop outlier cell lines
dat <- dat[, !colnames(dat) %in% exclude_cells]

# convert ensembl gene ids -> gene symbols
dat$symbol <- grch38$symbol[match(dat$gene_id, grch38$ensgene)]

#table(is.na(dat$symbol))
# 
# FALSE  TRUE 
# 24559   388 

# remove unmapped genes
dat <- dat[!is.na(dat$symbol), ]

#head(sort(table(dat$symbol), TRUE))
# KIR3DL3 KIR3DL2 KIR2DL4 KIR2DL1 KIR2DS4 KIR3DL1 
#      32      31      30      25      22      22 

#table(complete.cases(dat))
# FALSE  TRUE 
#   262 24297 

# remove genes with missing values
dat <- dat[complete.cases(dat), ]

# average multi-mapped gene symbols
dat <- dat %>%
  select(-gene_id) %>%
  group_by(symbol) %>%
  summarise_all(funs(mean))

# there are no zero-variance genes, so no filtering needed

write_tsv(dat, '/data/projects/nih/p3-models/data/raw/features/cnv/cnv_gene_max_scores_grch38.tsv')

