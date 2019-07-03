#!/bin/env Rscript
library(PharmacoGx)
library(Biobase)
library(tidyverse)
options(stringsAsFactors = FALSE)

# target drugs (choosing ones that were also measured in HMCL and appeared to show
# differential response across cell lines in that dataset)
target_drugs <- c('embelin', 'shikonin')

#avail <- availablePSets()
#avail[, 1:3]
#                       PSet.Name Dataset.Type Available.Molecular.Profiles
# CCLE_2013             CCLE_2013  sensitivity                 rna/mutation
# CCLE                       CCLE  sensitivity      rna/rnaseq/mutation/cnv
# GDSC_2013             GDSC_2013  sensitivity                 rna/mutation
# GDSC                       GDSC  sensitivity rna/rna2/mutation/fusion/cnv
# GDSC1000               GDSC1000  sensitivity                          rna
# gCSI                       gCSI  sensitivity                   rnaseq/cnv
# FIMM                       FIMM  sensitivity                             
# CTRPv2                   CTRPv2  sensitivity                             
# CMAP                       CMAP perturbation                          rna
# L1000_compounds L1000_compounds perturbation                          rna
# L1000_genetic     L1000_genetic perturbation                          rna
#ctrp <- downloadPSet('CTRPv2', saveDir = file.path('/data/public/human/', 'ctrpv2')) 

output_dir <- '/data/public/human/gdsc'
gdsc <- downloadPSet('GDSC', saveDir = output_dir) 

#
# rnaseq
#

# Extract RNA expression data to a matrix 
rna_eset <- summarizeMolecularProfiles(gdsc, mDataType = "rna") 

# collapse multi-mapped entries (120 / 11893)
rna_dat <- bind_cols(symbol = fData(rna_eset)$Symbol, as.data.frame(exprs(rna_eset)))

rna_dat <- rna_dat %>% 
  group_by(symbol) %>%
  summarize_all(mean) %>%
  ungroup

rna_dat <- rna_dat[!is.na(rna_dat$symbol), ]

#
# mutation data
#
mut_eset <- summarizeMolecularProfiles(gdsc, mDataType = "mutation", summary.stat = 'or') 
mut_dat <- as.data.frame(cbind(symbol = fData(mut_eset)$gene_name, exprs(mut_eset)))

# convert to numeric
for (cname in colnames(mut_dat)[-1]) {
  mut_dat[, cname] <- as.numeric(mut_dat[, cname])
}

#
# cnv data 
#
cnv_eset <- summarizeMolecularProfiles(gdsc, mDataType = "cnv")
cnv_dat <- as.data.frame(cbind(symbol = fData(cnv_eset)$Symbol, exprs(cnv_eset)))

# convert to numeric
for (cname in colnames(cnv_dat)[-1]) {
  cnv_dat[, cname] <- as.numeric(cnv_dat[, cname])
}

#
# drug response
#
#sensitivityMeasures(gdsc)
# [1] "ic50_published"  "auc_published"   "auc_recomputed"  "ic50_recomputed"

#sensitivityInfo(gdsc) %>% head
#                  cellid     drugid  drug.name nbr.conc.tested min.Dose.uM max.Dose.uM duration_h
# drugid_1_MC-CAR  MC-CAR  Erlotinib  ERLOTINIB               9 0.007812500      2.0000         72
# drugid_3_MC-CAR  MC-CAR  Rapamycin  RAPAMYCIN               9 0.000390625      0.1000         72
# drugid_5_MC-CAR  MC-CAR  Sunitinib  SUNITINIB               9 0.031250000      8.0000         72
# drugid_6_MC-CAR  MC-CAR PHA-665752  PHA665752               9 0.007812500      2.0000         72
# drugid_9_MC-CAR  MC-CAR     MG-132      MG132               9 0.003906250      1.0000         72
# drugid_11_MC-CAR MC-CAR paclitaxel PACLITAXEL               9 0.000400000      0.1024         72

#drug_auc <- summarizeSensitivityProfiles(gdsc, sensitivity.measure = 'auc_recomputed', 
#                                         drugs = target_drugs)

# save drug output
#for (drug_id in rownames(drug_auc)) {
#  auc_dat <- data.frame(cell_line = colnames(drug_auc), auc = drug_auc[drug_id, ]) 
#  write_tsv(auc_dat, file.path(output_dir, sprintf('%s_auc.tsv.gz', drug_id)))
#}

drug_ic50 <- summarizeSensitivityProfiles(gdsc, sensitivity.measure = 'ic50_recomputed', 
                                          drugs = target_drugs)

# find samples covered by all data types
#dsets <- list(rna_dat, mut_dat, cnv_dat)
#shared_sample_ids <- Reduce(intersect, lapply(dsets, function(x) { colnames(x)[-1] }))

# each data type covers all samples
# length(shared_sample_ids) == ncol(rna_dat) - 1
# [1] TRUE

#all(colnames(rna_dat) == colnames(mut_dat)) 
# [1] TRUE

#all(colnames(rna_dat) == colnames(cnv_dat)) 
# [1] TRUE

for (drug_id in rownames(drug_ic50)) {
  ic50_dat <- data.frame(cell_line = colnames(drug_ic50), ic50 = drug_ic50[drug_id, ]) 

  # construction a 1 x <num samples> dataframe with response data
  ic50_dat <- t(as.data.frame(setNames(ic50_dat$ic50, ic50_dat$cell_line))) 
  ic50_dat <- cbind(data.frame(response='ic50_recomputed'), ic50_dat)

  # all(colnames(ic50_dat)[-1] == colnames(rna_dat)[-1])
  # [1] TRUE

  write_tsv(ic50_dat, file.path(output_dir, sprintf('GDSC_%s_ic50_recomputed.tsv.gz', drug_id)))
}

# save sample metadata
write_tsv(pData(rna_eset), file.path(output_dir, 'GDSC_samples.tsv'))

# save feature data
write_tsv(mut_dat, file.path(output_dir, 'GDSC_mut.tsv.gz'))
write_tsv(rna_dat, file.path(output_dir, 'GDSC_rna.tsv.gz'))
write_tsv(cnv_dat, file.path(output_dir, 'GDSC_cnv.tsv.gz'))

