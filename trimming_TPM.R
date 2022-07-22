


##### using EdgeR to filter genes with low expression #####


rm(list = ls())
library('tidyverse')
library("ggpubr")
library('edgeR')

direc <- 'C:\\Users\\Ben\\OneDrive - Duke University\\bilbo_lab\\alexis_sequencing\\'

## re the edgeR users guide we should calculate CPM on the raw counts and then trim
## https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


### reading in and cleaning raw counts ###
count_df <- read.csv(paste(direc, 'counts\\count_mtx_FINAL_w_manual_exclusions.csv', sep = '\\'))
######## NEED TO PULL IN THE COUNT MATRIX WITH SAMPLES REMOVED ########
count_df <- count_df[!duplicated(count_df$hgnc_symbol), ]
rownames(count_df) <- count_df$hgnc_symbol


### separate brain and placenta samples ###

placenta <- count_df[, grep('placenta', colnames(count_df))]
brain <- count_df[, grep('brain', colnames(count_df))]


## analyzing placenta
placenta_cts_cpm <- cpm(placenta)


placenta_cts_cpm_dge <- DGEList(counts=placenta_cts_cpm)
keep <- filterByExpr(placenta_cts_cpm_dge, min.count = 1)
placenta_cts_cpm_dge <- placenta_cts_cpm_dge[keep, keep.lib.sizes=FALSE]
plac_output <- placenta_cts_cpm_dge$counts

write.csv(rownames(plac_output), paste(direc, 'counts\\placenta_genes_to_keep.csv', sep = '\\'))

placenta$gene <- row.names(placenta)
plac_output <- data.frame(plac_output)
plac_output$gene <- row.names(plac_output)
write.csv(merge(placenta, plac_output, by.x = 'gene', by.y = 'gene'), paste(direc, 'counts\\trimmed_placenta_counts.csv', sep = '\\'))



## analyzing brain
brain_cts_cpm <- cpm(brain)
brain_cts_cpm_dge <- DGEList(counts=brain_cts_cpm)
keep <- filterByExpr(brain_cts_cpm_dge, min.count = 1)
brain_cts_cpm_dge <- brain_cts_cpm_dge[keep, keep.lib.sizes=FALSE]
brain_output <- brain_cts_cpm_dge$counts

write.csv(rownames(brain_output), paste(direc, 'counts\\brain_genes_to_keep.csv', sep = '\\'))


brain$gene <- row.names(brain)
brain_output <- data.frame(brain_output)
brain_output$gene <- row.names(brain_output)
write.csv(merge(brain, brain_output, by.x = 'gene', by.y = 'gene'), paste(direc, 'counts\\trimmed_brain_counts.csv', sep = '\\'))




######## TRYING THE EDGER TRIMMING WITH TPM INSTEAD OF CPM ########
brain_tpm <- read.csv(paste(direc, 'counts\\brain_tpm_NO_trim.csv', sep = '\\'))
brain_tpm_dge <- DGEList(counts=brain_tpm[5:38])
row.names(brain_tpm_dge) <- brain_tpm$hgnc_symbol
keep <- filterByExpr(brain_tpm_dge, min.count = 1)
brain_tpm_dge <- brain_tpm_dge[keep, keep.lib.sizes=FALSE]
brain_tpm_dge_output <- brain_tpm_dge$counts

write.csv(brain_tpm_dge_output, paste(direc, 'counts\\brain_tpm_TRIM.csv', sep = '\\'))

placenta_tpm <- read.csv(paste(direc, 'counts\\placenta_tpm_NO_trim.csv', sep = '\\'))
placenta_tpm_dge <- DGEList(counts=placenta_tpm[5:38])
row.names(placenta_tpm_dge) <- placenta_tpm$hgnc_symbol
keep <- filterByExpr(placenta_tpm_dge, min.count = 1)
placenta_tpm_dge <- placenta_tpm_dge[keep, keep.lib.sizes=FALSE]
placenta_tpm_dge_output <- placenta_tpm_dge$counts

write.csv(placenta_tpm_dge_output, paste(direc, 'counts\\placenta_tpm_TRIM.csv', sep = '\\'))





