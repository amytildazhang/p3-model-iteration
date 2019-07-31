#!/bin/env Rscript
library(tidyverse)
options(stringsAsFactors = FALSE)

data_dir <- '/data/projects/nih/p3/drug_response/curves/1.2' 
output_dir <- '/data/projects/nih/p3-models/data/raw/response'

target_drugs <- c('A-443654', 'Adenosine', 'Aminopterin', 'CGS-21680', 
                  'Cyproheptadine hydrochloride', 'Embelin', 'EMD-1214063', 'Finasteride',
                  'Flutamide', 'GM-6001', 'Mubritinib', 'PD-173955 Analogue 1',
                  'RGB-286147', 'SB-505124', 'Secoisolariciresinol', 'Shikonin',
                  'Vinflunine', 'WHI-P97', 'WZ-4002')

exclude_cells <- c('KMS21BM_JCRB', 'Karpas417_ECACC', 'OCIMY1_PLB', 'OPM2_DSMZ')

# iterate over different versions of drug response data
for (infile in Sys.glob(file.path(data_dir, '*.csv'))) {
  # load data
  dat <- read_csv(infile) %>%
    filter(drug_name %in% target_drugs) %>%
    filter(!cell_line %in% exclude_cells) %>%
    select(cell_line, drug_id, ac50, lac50, DATA7, DATA8)

  # base output filename
  file_prefix <- tools::file_path_sans_ext(basename(infile))

  # iterate over response types, and for each drug / response, save to a separate file
  for (response in c('ac50', 'lac50', 'DATA7', 'DATA8')) {
    
    dat %>%
      select(drug_id, cell_line, response) %>%
      group_by(drug_id) %>%
      group_walk(~write_tsv(.x, file.path(output_dir, sprintf("%s_%s_%s.tsv", .y$drug_id, file_prefix, response))))
  }
}

