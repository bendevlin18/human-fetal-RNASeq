


library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)
library(plotly)
library(ggpubr)
library(DT)
library(tidyverse)

setwd('C:\\Users\\Ben\\OneDrive - Duke University\\bilbo_lab\\alexis_sequencing\\shiny_app')


###### TISSUE TPM GENE MATRICES ######

brain_tpm_df <- read.csv('brain_tpm_w_exclusions.csv')
placenta_tpm_df <- read.csv('placenta_tpm_w_exclusions.csv')

brain_genes <- brain_tpm_df$hgnc_symbol

z <- t(brain_tpm_df)
## batch through and clean colnames
for (cols in colnames(brain_tpm_df)) {
  
  if (grepl('tpm', cols, fixed=TRUE)) {
    
    colnames(brain_tpm_df)[colnames(brain_tpm_df) == cols] <- str_split(str_split(cols, 'X')[[1]][2], '_tpm')[[1]][1]
    
  }
}


###### TISSUE AND SAMPLE METADATA #######

metadata <- read.csv('human_seq_metadata.csv', fileEncoding = 'UTF-8-BOM')


##################################################

metadata$sampleID <- metadata$ï..Sample.ID


for (sample in metadata$sampleID) {
  
  metadata$sampleID[metadata$sampleID == sample] <- str_split(sample, "-")[[1]][2]
  
}

metadata <- metadata[!grepl('y', metadata$exclude),]

#SLC6A4, TLR4, TPH1, TH2, ido1, ido2, MAOA, RPS4Y1, nlgn4y, CD180, ccl2, cd14, cd163, TPH2, tmem119, HTR1A, HTR2A, DDC, emr1 (ADGRE1), RRM2, cd68

## change this variable to explore different genes
## make sure to use the common name, in ALL CAPS
## run from here down
gene_of_interest <- 'TLR5'




## creating the gene of interest (goi) dataframe
goi_df <- data.frame(brain_tpm_df[brain_tpm_df$hgnc_symbol == gene_of_interest, ])
goi_df <- data.frame(t(goi_df))
goi_df$sampleIDs <- rownames(goi_df)
rownames(goi_df) <- NULL
goi_df <- goi_df[4:81,]
colnames(goi_df)[1] <- gene_of_interest

## looping to extract the sampleID column and the region column
goi_df$sampleID <- goi_df$sampleIDs
goi_df$region <- goi_df$sampleIDs
for (sample in goi_df$sampleIDs) {
  
  
  
  goi_df$sampleID[goi_df$sampleIDs == sample] <- str_split(str_split(sample, 'X')[[1]][2], '_')[[1]][1]
  
  goi_df$region[goi_df$sampleIDs == sample] <- str_split(str_split(sample, 'X')[[1]][2], '_')[[1]][2]
  
}


## merging goi_df with the metadata
merged <- merge(goi_df, metadata, by.x="sampleID", by.y="sampleID")
merged <- transform(merged, `TLR5` = as.numeric(`TLR5`))
merged$sex <- factor(merged$sex, levels = c('Male', 'Female'))


ggplotly(ggscatter(merged, x = 'trig', y = 'TLR5', add = 'reg.line', conf.int = TRUE, cor.method = "pearson")+
  ggtitle(paste(gene_of_interest, 'Brain', sep = ' '))+
  facet_wrap(~sex)+
  xlab('Log Triglycerides (mg/dL)')+
  ylab(paste(gene_of_interest, '(TPM)', sep = ' '))+
  theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5)))

