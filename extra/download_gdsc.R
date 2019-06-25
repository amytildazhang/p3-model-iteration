#!/bin/env Rscript
library(PharmacoGx)
library(tidyverse)
options(stringsAsFactors = FALSE)

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
gdsc_eset <- summarizeMolecularProfiles(gdsc, mDataType = "rna") 

# save sample metadata
write_tsv(pData(gdsc_eset), file.path(output_dir, 'GDSC_rna_samples.tsv'))

# collapse multi-mapped entries (120 / 11893)
expr_dat <- bind_cols(symbol = fData(gdsc_eset)$Symbol, as.data.frame(exprs(gdsc_eset)))

expr_dat <- expr_dat %>% 
  group_by(symbol) %>%
  summarize_all(mean) %>%
  ungroup

expr_dat <- expr_dat[!is.na(expr_dat$symbol), ]
write_tsv(expr_dat, file.path(output_dir, 'GDSC_rna_expr.tsv.gz'))

#
# mutation data
#
gdsc_eset <- summarizeMolecularProfiles(gdsc, mDataType = "mutation", summary.stat = 'or') 

mut <- as.data.frame(exprs(gdsc_eset))

for (cname in colnames(mut)) {
  mut[, cname] <- as.numeric(mut[, cname])
}

write_tsv(mut, file.path(output_dir, 'GDSC_mutations.tsv.gz'))

#
# cnv data 
#
gdsc_eset <- summarizeMolecularProfiles(gdsc, mDataType = "cnv")

cnv <- as.data.frame(exprs(gdsc_eset))

write_tsv(cnv, file.path(output_dir, 'GDSC_cnv.tsv.gz'))

#
# drug response
#
#sensitivityMeasures(gdsc)
# [1] "ic50_published"  "auc_published"   "auc_recomputed"  "ic50_recomputed"

sensitivityInfo(gdsc) %>% head
#                  cellid     drugid  drug.name nbr.conc.tested min.Dose.uM max.Dose.uM duration_h
# drugid_1_MC-CAR  MC-CAR  Erlotinib  ERLOTINIB               9 0.007812500      2.0000         72
# drugid_3_MC-CAR  MC-CAR  Rapamycin  RAPAMYCIN               9 0.000390625      0.1000         72
# drugid_5_MC-CAR  MC-CAR  Sunitinib  SUNITINIB               9 0.031250000      8.0000         72
# drugid_6_MC-CAR  MC-CAR PHA-665752  PHA665752               9 0.007812500      2.0000         72
# drugid_9_MC-CAR  MC-CAR     MG-132      MG132               9 0.003906250      1.0000         72
# drugid_11_MC-CAR MC-CAR paclitaxel PACLITAXEL               9 0.000400000      0.1024         72

# target drugs (choosing ones that were also measured in HMCL and appeared to show
# differential response across cell lines in that dataset)
target_drugs <- c('embelin', 'shikonin')

drug_auc <- summarizeSensitivityProfiles(gdsc, sensitivity.measure = 'auc_recomputed', 
                                         drugs = target_drugs)

drug_ic50 <- summarizeSensitivityProfiles(gdsc, sensitivity.measure = 'ic50_recomputed', 
                                         drugs = target_drugs)

# save drug output
for (drug_id in rownames(drug_auc)) {
  auc_dat <- data.frame(cell_line = colnames(drug_auc), auc = drug_auc[drug_id, ]) 
  write_tsv(auc_dat, file.path(output_dir, sprintf('%s_auc.tsv.gz', drug_id)))

  ic50_dat <- data.frame(cell_line = colnames(drug_ic50), ic50 = drug_ic50[drug_id, ]) 
  write_tsv(ic50_dat, file.path(output_dir, sprintf('%s_ic50.tsv.gz', drug_id)))
}


