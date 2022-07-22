

rm(list = ls())
library('biomaRt')
library('EDASeq')
library('tidyverse')

## post that outlined this workflow
## https://stackoverflow.com/questions/28543517/how-can-i-convert-ensembl-id-to-gene-symbol-in-r

## grabbing homo sapien gene ensembl ID matches on BiomaRt
human <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

## reading in the counts file
###

## list all the files in given directory
count_direc <- 'C:\\Users\\Ben\\OneDrive - Duke University\\bilbo_lab\\alexis_sequencing\\counts'
files <- list.files(count_direc)

## load in first file for a template df
df <- read.table(paste(count_direc, files[1], sep = '\\'), fill = TRUE, header = 1)
df <- select(df, c('Length', 'Geneid'))


## read in all the files in the directory and adding them to the template df
for (file in files) {
  if (grepl('.txt', file, fixed=TRUE)) {
    df[file] <- read.table(paste(count_direc, file, sep = '\\'), fill = TRUE, header = 1)[7]
    print(paste('Reading..', file))
  }
}

## batching through the column names to get rid of the "counts.txt" suffix
for (cols in colnames(df)) {
  
  if (grepl('.txt', cols, fixed=TRUE)) {
    colnames(df)[colnames(df) == cols] <- c(str_split(c(cols), '_counts.txt')[[1]][1])
    
  }
}

# 
###

## creating gene list with hgnc symbols from BiomaRt
gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = df$Geneid, mart = human)

## merging hgnc symbols with count table using ensembl IDs
merged <- merge(df, gene_list, by.x="Geneid", by.y="ensembl_gene_id")


########## ALL OF THIS IS DEPRECATED ##########
## getting the gene lengths using EDASeq and biomart database
## this took ~45mins, so I saved it to a file to read in later
##ensembl_list <- c(merged$Geneid)
#gene_length <- getGeneLengthAndGCContent(ensembl_list, "hsa")
#write.csv(gene_length, 'C:\\Users\\Ben\\Documents\\alexis_sequencing\\metadata\\gene_lengths.csv')
##gene_length <- read.csv('C:\\Users\\Ben\\Documents\\alexis_sequencing\\metadata\\gene_lengths.csv')

## merging together the gene lengths, common gene names, and all counts for all samples
##final_df <- merge(merged, gene_length, by.x = 'Geneid', by.y = 'X')
##final_df$gc <- NULL

write.csv(merged, paste(count_direc, 'count_mtx_FINAL.csv', sep = '\\'))



####################################################
####################################################
############## Now time to create TPM ##############
####################################################
####################################################



## TPM is very similar to RPKM and FPKM. The only difference is the order of operations. Here’s how you calculate TPM:
## eq from https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/ 

count_df <- read.csv(paste(count_direc, 'count_mtx_FINAL_w_manual_exclusions.csv', sep = '\\'))
count_df$X <- NULL
count_df <- count_df[!duplicated(count_df$hgnc_symbol), ]


### separate brain and placenta samples ###


#### we're first going to do the operation on all the placenta samples ####

placenta <- count_df[, grep('placenta', colnames(count_df))]
placenta$hgnc_symbol <- count_df$hgnc_symbol
placenta$Length <- count_df$Length
placenta$Geneid <- count_df$Geneid

## batching through the column names to calculate TPM for each sample
for (cols in colnames(placenta)) {
  
  if (grepl('X', cols, fixed=TRUE)) {
    
    print(cols)
    #Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
    lengths <- placenta$Length
    rpk <- placenta[cols] / lengths
    
    #Count up all the RPK values in a sample and divide this number by 1,000,000. This is your per million scaling factor.
    per_mil <- sum(rpk, na.rm=T) / 1000000
    
    #Divide the RPK values by the per million scaling factor. This gives you TPM.
    tpm <- rpk/per_mil
    
    
    placenta[paste(cols, 'tpm', sep = '_')] <- tpm

    
  }
}

## now i wanna batch through and get rid of the counts
placenta_tpm_df <- dplyr::select(placenta, c('Length', 'Geneid', 'hgnc_symbol'))
for (cols in colnames(placenta)) {
  
  if (grepl('tpm', cols, fixed=TRUE)) {
    
    print(cols)
    placenta_tpm_df[cols] <- placenta[cols]
    
    
  }
}

write.csv(placenta_tpm_df, paste(count_direc, 'placenta_tpm_NO_trim.csv', sep = '\\'))
placenta_genes <- read.csv(paste(count_direc, 'placenta_genes_to_keep.csv', sep = '\\'))
placenta_tpm_final <- merge(placenta_tpm_df, placenta_genes, by.x = 'hgnc_symbol', by.y = 'x')
write.csv(placenta_tpm_final, paste(count_direc, 'placenta_tpm_w_exclusions.csv', sep = '\\'))




#### now we are going to do the operation on all the brain samples ####

brain <- count_df[, grep('brain', colnames(count_df))]
brain$hgnc_symbol <- count_df$hgnc_symbol
brain$Length <- count_df$Length
brain$Geneid <- count_df$Geneid

## batching through the column names to calculate TPM for each sample
for (cols in colnames(brain)) {
  
  if (grepl('X', cols, fixed=TRUE)) {
    
    print(cols)
    #Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
    lengths <- brain$Length
    rpk <- brain[cols] / lengths
    
    #Count up all the RPK values in a sample and divide this number by 1,000,000. This is your per million scaling factor.
    per_mil <- sum(rpk, na.rm=T) / 1000000
    
    #Divide the RPK values by the per million scaling factor. This gives you TPM.
    tpm <- rpk/per_mil
    
    
    brain[paste(cols, 'tpm', sep = '_')] <- tpm
    
    
  }
}

## now i wanna batch through and get rid of the counts
brain_tpm_df <- dplyr::select(placenta, c('Length', 'Geneid', 'hgnc_symbol'))
for (cols in colnames(brain)) {
  
  if (grepl('tpm', cols, fixed=TRUE)) {
    
    print(cols)
    brain_tpm_df[cols] <- brain[cols]
    
    
  }
}

write.csv(brain_tpm_df, paste(count_direc, 'brain_tpm_NO_trim.csv', sep = '\\'))
brain_genes <- read.csv(paste(count_direc, 'brain_genes_to_keep.csv', sep = '\\'))
brain_tpm_final <- merge(brain_tpm_df, brain_genes, by.x = 'hgnc_symbol', by.y = 'x')
write.csv(brain_tpm_final, paste(count_direc, 'brain_tpm_w_exclusions.csv', sep = '\\'))






