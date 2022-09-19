

######################################
###                                ###
###     global.R file for          ###
###       human-fetal-seq          ###
###        application             ###
###                                ###
######################################



##### Load in necessary packages #####

library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)
library(plotly)
library(ggpubr)
library(DT)
library(tidyverse)
library(shinycssloaders)
library(shinythemes)



###### TISSUE TPM GENE MATRICES ######
brain_tpm_df <- read.csv('brain_tpm_TRIM_MBX.csv')
placenta_tpm_df <- read.csv('placenta_tpm_TRIM.csv')

brain_genes <- brain_tpm_df$hgnc_symbol
placenta_genes <- placenta_tpm_df$hgnc_symbol
choice_genes <- unique(list(brain_tpm_df$hgnc_symbol, placenta_tpm_df$hgnc_symbol))

##z <- t(brain_df)
## batch through and clean colnames for the brain
##for (cols in colnames(brain_df)) {
  
##  if (grepl('tpm', cols, fixed=TRUE)) {
    
##    colnames(brain_df)[colnames(brain_df) == cols] <- str_split(str_split(cols, 'X')[[1]][2], '_tpm')[[1]][1]
    
##
##}


##z <- t(placenta_df)
## batch through and clean colnames for the placenta
##for (cols in colnames(placenta_df)) {
  
##  if (grepl('tpm', cols, fixed=TRUE)) {
    
##    colnames(placenta_df)[colnames(placenta_df) == cols] <- str_split(str_split(cols, 'X')[[1]][2], '_tpm')[[1]][1]
    
##  }
##}



###### TISSUE AND SAMPLE METADATA #######

metadata <- read.csv('human_seq_metadata.csv', fileEncoding = 'UTF-8-BOM')

### prepping metadata
metadata$sampleID <- metadata$sampleID
for (sample in metadata$sampleID) {
  metadata$sampleID[metadata$sampleID == sample] <- str_split(sample, "-")[[1]][2]
}
metadata <- metadata[!grepl('y', metadata$exclude),]



corr_plot_settings <- list(
  xlab('Log Maternal Triglycerides (mg/dL)'),
  ylab('TPM (Transcripts per million)'),
  theme_classic(),
  theme(rect = element_rect(fill = 'transparent'), 
        text = element_text(size = 13),
        plot.title = element_text(size = 16, hjust = 0.5, face = "italic"),
        plot.margin = margin(t=2, r=1, b=0, l=1, unit="cm"),
        axis.line = element_line(size = .75),
        axis.ticks = element_line(size = .6),
        axis.text.x = element_text(hjust = 1),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color = NA, col = 0),
        axis.title.y = element_text(margin = margin(t = .05, r = 0.05, b = 0.05, l = 0.05), size=10),
        axis.title.x = element_text(margin = margin(t = .05, r = 0.05, b = 0.05, l = 0.05), size = 10),
        strip.background = element_rect(color="transparent", fill="transparent", linetype="solid")))

##next added for Sex Diffs plots

sem <- function(x) sqrt(var(x)/length(x))

box_plot_theme <- theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(fill = "transparent"), 
                        plot.background = element_rect(fill = "transparent"),
                        panel.border = element_rect(color="black", fill="NA", size=0.75),
                        legend.background = element_rect(fill = "transparent"),
                        plot.title = element_text(hjust = 0.5, face = "italic", size=11),
                        plot.margin = margin(t=.5, r=1, b=0, l=1, unit="cm"))

yy <- list(title = "TPM")


