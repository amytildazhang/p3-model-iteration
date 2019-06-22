#!/bin/env Rscript
gene_set_files <- c(bader_files, dsigdb_files, enrichr_files, gdsc_files, msigdb_files)

gene_sets <- lapply(gene_set_files, function(x) { geneIds(getGmt(x)) })

names(gene_sets) <- tools::file_path_sans_ext(basename(gene_set_files)) 
