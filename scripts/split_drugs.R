#!/bin/env Rscript
#
# Splits drug dataset csv into per-drug files for target drugs
#
library(tidyverse)

# load drug dataset
dat <- read_tsv(snakemake@input[[1]]) 

# extract target drugs and generate corresponding outputs
for (i in 1:length(snakemake@output)) {
  drug <- snakemake@config$target_drugs[i]

  dat %>%
    filter(drug_name == drug) %>%
    filter(!cell_line %in% snakemake@config$exclude_cell_lines) %>%
    write_tsv(snakemake@output[[i]])
}
