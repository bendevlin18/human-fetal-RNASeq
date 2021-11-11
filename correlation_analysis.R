

##### Plotting correlations between triglyceride and/or serotonin levels and genes of interest #####


rm(list = ls())
library('tidyverse')
library("ggpubr")

z <- theme(rect = element_rect(fill = 'transparent'), 
           text = element_text(size = 35, family = 'sans'), 
           axis.line = element_line(size = 2.5),
           axis.ticks = element_line(size = 1.7),
           axis.ticks.length = unit(.5, 'cm'),
           panel.background = element_rect(fill = "transparent"), # bg of the panel
           plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
           panel.grid.major = element_blank(), # get rid of major grid
           panel.grid.minor = element_blank(), # get rid of minor grid
           legend.background = element_rect(color = NA, col = 0),  # get rid of legend bg and outline
           legend.title = element_blank(),
           axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 30),
           axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
           strip.background = element_rect(color="transparent", fill="transparent", size=5, linetype="solid"),
           panel.spacing = unit(.5, 'cm'),
           axis.text.x = element_text(angle = 45, hjust = 1, size = 35),
           axis.text.y = element_text(size = 2),
           plot.title = element_text(hjust = 0.5, size = 40),
           legend.key.size = unit(1, 'cm'))

direc <- 'C:\\Users\\Ben\\OneDrive - Duke University\\bilbo_lab\\alexis_sequencing\\'


### reading in and cleaning placenta TPM dataframe ###
tpm_df <- read.csv(paste(direc, 'counts\\placenta_tpm_w_exclusions.csv', sep = '\\'))

z <- t(tpm_df)
## batch through and clean colnames
for (cols in colnames(tpm_df)) {
  
  if (grepl('tpm', cols, fixed=TRUE)) {

    colnames(tpm_df)[colnames(tpm_df) == cols] <- str_split(str_split(cols, 'X')[[1]][2], '_tpm')[[1]][1]
    
  }
}



### reading in and cleaning metadata dataframe ###
metadata <- read.csv(paste(direc, 'metadata/human_seq_metadata.csv', sep = '/'))


metadata$sampleID <- metadata$ï..Sample.ID
for (sample in metadata$sampleID) {
  
  metadata$sampleID[metadata$sampleID == sample] <- str_split(sample, "-")[[1]][2]
  
}

metadata <- metadata[!grepl('y', metadata$exclude),]


#SLC6A4, TLR4, TPH1, TH2, ido1, ido2, MAOA, RPS4Y1, nlgn4y, CD180, ccl2, cd14, cd163, CSF1Rr, TPH2, tmem119, HTR1A, HTR2A, DDC, emr1 (ADGRE1), APOE, cd68

## change this variable to explore different genes
## make sure to use the common name, in ALL CAPS
gene_of_interest <- 'CSF1R'




## creating the gene of interest (goi) dataframe
goi_df <- data.frame(tpm_df[tpm_df$hgnc_symbol == gene_of_interest, ])
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

merged <- transform(merged, `CSF1R` = as.numeric(`CSF1R`))
merged$sex <- factor(merged$sex, levels = c('Male', 'Female'))


## now plot placenta
merged_plac <- subset(merged, merged$region == 'placenta')
ggscatter(merged_plac, x = 'trig', y = 'CSF1R',
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson"
          )+
  facet_wrap(~sex)+
  ggtitle(paste(gene_of_interest, 'Placenta', sep = ' '))+
  xlab('Log Triglycerides (mg/dL)')+
  ylab(paste(gene_of_interest, '(TPM)', sep = ' '))+
  theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5))


ggsave(paste(gene_of_interest, 'plac_corr_plot.svg', sep = '_'))


#################################################################################################################

### reading in and cleaning placenta TPM dataframe ###
tpm_df <- read.csv(paste(direc, 'counts\\brain_tpm_w_exclusions.csv', sep = '\\'))

z <- t(tpm_df)
## batch through and clean colnames
for (cols in colnames(tpm_df)) {
  
  if (grepl('tpm', cols, fixed=TRUE)) {
    
    colnames(tpm_df)[colnames(tpm_df) == cols] <- str_split(str_split(cols, 'X')[[1]][2], '_tpm')[[1]][1]
    
  }
}



### reading in and cleaning metadata dataframe ###
metadata <- read.csv(paste(direc, 'metadata/human_seq_metadata.csv', sep = '/'))


metadata$sampleID <- metadata$ï..Sample.ID
for (sample in metadata$sampleID) {
  
  metadata$sampleID[metadata$sampleID == sample] <- str_split(sample, "-")[[1]][2]
  
}

metadata <- metadata[!grepl('y', metadata$exclude),]


#SLC6A4, TLR4, TPH1, TH2, ido1, ido2, MAOA, RPS4Y1, nlgn4y, CD180, ccl2, cd14, cd163, CSF1Rr, TPH2, tmem119, HTR1A, HTR2A, DDC, emr1 (ADGRE1), APOE, cd68

## change this variable to explore different genes
## make sure to use the common name, in ALL CAPS
gene_of_interest <- 'CSF1R'




## creating the gene of interest (goi) dataframe
goi_df <- data.frame(tpm_df[tpm_df$hgnc_symbol == gene_of_interest, ])
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

merged <- transform(merged, `CSF1R` = as.numeric(`CSF1R`))
merged$sex <- factor(merged$sex, levels = c('Male', 'Female'))


## plot brain first
merged_brain <- subset(merged, merged$region == 'brain')
ggscatter(merged_brain, x = 'trig', y = 'CSF1R',
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson"
)+
  facet_wrap(~sex)+
  ggtitle(paste(gene_of_interest, 'Brain', sep = ' '))+
  xlab('Log Triglycerides (mg/dL)')+
  ylab(paste(gene_of_interest, '(TPM)', sep = ' '))+
  theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5))


ggsave(paste(gene_of_interest, 'brain_corr_plot.svg', sep = '_'))



