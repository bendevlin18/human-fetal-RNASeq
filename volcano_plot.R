

rm(list = ls())
library(tidyverse)
library(ggpubr)
library(reshape2)
library(ggbeeswarm)
library(DESeq2)
library(openxlsx)
library(EnhancedVolcano)
library(ggbreak)
library(mltools)
library(ggrepel)

setwd('C:\\Users\\Ben\\OneDrive - Duke University\\bilbo_lab\\alexis_sequencing\\analysis\\pval_cutoffs\\w_TPM_trimming')
##setwd('/Users/bendevlin/OneDrive - Duke University/bilbo_lab/alexis_sequencing/analysis')
m_brain <- read.csv('pearson_corr_w_triglycerides_m_brain.csv')

m_brain$log10_pval <- -log10(m_brain$pval)

### binning the correlation coefficients into 100 bins
bins <- bin_data(m_brain$CorrCoef, bins = 50, returnDT = TRUE)

### creating a new dataframe out of these bins which assigns a numeric value to each bin
binz <- unique(bins$Bin)
final_bins <- data.frame(binz)
final_bins$index <- seq(1, 50, by=1)

bins_final <- merge(bins, final_bins, by.x = 'Bin', by.y = 'binz')
### merging the bin dataframe with the final dataframe
new_df <- merge(m_brain, bins_final, by.x = 'CorrCoef', by.y = 'BinVal')


avg <- new_df %>%
  group_by(index) %>%
  summarize(binned_corr_coef = mean(CorrCoef))

final_merged <- merge(new_df, avg, by = 'index')

ggplot(data = final_merged, aes(x = binned_corr_coef, y =log10_pval, label = Gene))+
  geom_quasirandom(data = subset(final_merged, log10_pval < -log10(.05)), size = 1, color = 'black', width = .2, dodge.width = .5)+
  geom_hline(yintercept = -log10(.05), linetype = 'dotted', size = 1.5)+
  ## this is the line where we can add the genes we want to highlight
  geom_label_repel(data = subset(final_merged, Gene == 'TEKT5' | Gene == 'MIR93'), box.padding = 1, max.overlaps = Inf)+
  geom_quasirandom(data = subset(final_merged, log10_pval > -log10(.05)), size = 2, color = 'red', width = .2, dodge.width = .5)+
  #geom_label_repel(data = subset(final_merged, log10_pval > 3 & CorrCoef > .4 | CorrCoef < -.5),
  #box.padding   = 0.35, 
  #point.padding = 0.5,
  #segment.color = 'grey50')+
  ggtitle('Male Brain')+
  xlab('Binned CorrCoef')+
  ylab('Log10_Pval')+
  theme_bw()+
  scale_y_continuous(limits = c(0,6), expand = c(0, 0))+
  theme(plot.title = element_text(hjust = 0.5))

write.csv(final_merged, 'male_brain_w_bins.csv')

ggsave('male_brain_FINAL_volcano_plot_for_Alexis_paper.svg')



#### reading in male placenta

m_brain <- read.csv('pearson_corr_w_triglycerides_m_plac.csv')

m_brain$log10_pval <- -log10(m_brain$pval)


### binning the correlation coefficients into 100 bins
bins <- bin_data(m_brain$CorrCoef, bins = 50, returnDT = TRUE)

### creating a new dataframe out of these bins which assigns a numeric value to each bin
binz <- unique(bins$Bin)
final_bins <- data.frame(binz)
final_bins$index <- seq(1, 50, by=1)


bins_final <- merge(bins, final_bins, by.x = 'Bin', by.y = 'binz')
### merging the bin dataframe with the final dataframe
new_df <- merge(m_brain, bins_final, by.x = 'CorrCoef', by.y = 'BinVal')
new_df <- new_df[!duplicated(new_df$Gene), ]

avg <- new_df %>%
  group_by(index) %>%
  summarize(binned_corr_coef = mean(CorrCoef))

final_merged <- merge(new_df, avg, by = 'index')

ggplot(data = final_merged, aes(x = binned_corr_coef, y =log10_pval, label = Gene))+
  geom_quasirandom(data = subset(final_merged, log10_pval > -log10(.05)), size = 2, color = 'red', width = .2, dodge.width = .5)+
  geom_quasirandom(data = subset(final_merged, log10_pval < -log10(.05)), size = 1, color = 'black', width = .2, dodge.width = .5)+
  geom_hline(yintercept = -log10(.05), linetype = 'dotted', size = 1.5)+
  geom_label_repel(data = subset(final_merged, Gene == 'NFE2' | Gene == 'TUBB1' | Gene == 'FAM210B' | Gene == 'TSPO2' | 
                                  Gene == 'CSF1' | Gene == 'CD163' | Gene == 'CD14' | 
                                  Gene == 'CCL2' | Gene == 'IL17RA' | Gene == 'CSF1R' | Gene == 'NFKB1'
                                ), box.padding = 2, max.overlaps = Inf)+
  #box.padding   = 0.35, 
  #point.padding = 0.5,
  #segment.color = 'grey50')+
  ggtitle('Male Placenta')+
  xlab('Binned CorrCoef')+
  ylab('Log10_Pval')+
  theme_bw()+
  scale_y_continuous(limits = c(0,6), expand = c(0, 0))+
  theme(plot.title = element_text(hjust = 0.5))

write.csv(final_merged, 'male_placenta_w_bins.csv')

ggsave('male_plac_FINAL_volcano_plot_for_Alexis_paper.svg')


#### reading in female placenta

m_brain <- read.csv('pearson_corr_w_triglycerides_f_plac.csv')

m_brain$log10_pval <- -log10(m_brain$pval)


### binning the correlation coefficients into 100 bins
bins <- bin_data(m_brain$CorrCoef, bins = 50, returnDT = TRUE)

### creating a new dataframe out of these bins which assigns a numeric value to each bin
binz <- unique(bins$Bin)
final_bins <- data.frame(binz)
final_bins$index <- seq(1, 50, by=1)


bins_final <- merge(bins, final_bins, by.x = 'Bin', by.y = 'binz')
### merging the bin dataframe with the final dataframe
new_df <- merge(m_brain, bins_final, by.x = 'CorrCoef', by.y = 'BinVal')
##drop the duplicates from this dataframe bc there is a weird error i cannot figure out
new_df <- new_df[!duplicated(new_df$Gene), ]

avg <- new_df %>%
  group_by(index) %>%
  summarize(binned_corr_coef = mean(CorrCoef))

final_merged <- merge(new_df, avg, by = 'index')

ggplot(data = final_merged, aes(x = binned_corr_coef, y =log10_pval, label = Gene))+
  geom_quasirandom(data = subset(final_merged, log10_pval > -log10(.05)), size = 2, color = 'red', width = .2, dodge.width = .5)+
  geom_quasirandom(data = subset(final_merged, log10_pval < -log10(.05)), size = 1, color = 'black', width = .2, dodge.width = .5)+
  geom_hline(yintercept = -log10(.05), linetype = 'dotted', size = 1.5)+
  geom_label_repel(data = subset(final_merged, Gene == 'PECAM1' | Gene == 'SPP1' | Gene == 'PDGFA' | Gene == 'VCAM' | 
                                   Gene == 'CD36' | Gene == 'CXCL1' | Gene == 'NFKBIA' | 
                                   Gene == 'TNFAIP3' | Gene == 'MCL1' | Gene == 'DUSP1'
  ), box.padding = 2, max.overlaps = Inf)+
  #geom_label_repel(data = subset(final_merged, log10_pval > 3 & CorrCoef > .4 | CorrCoef < -.5),
  #box.padding   = 0.35, 
  #point.padding = 0.5,
  #segment.color = 'grey50')+
  ggtitle('Female Placenta')+
  xlab('Binned CorrCoef')+
  ylab('Log10_Pval')+
  theme_bw()+
  scale_y_continuous(limits = c(0,6), expand = c(0, 0))+
  theme(plot.title = element_text(hjust = 0.5))


write.csv(final_merged, 'female_plac_w_bins.csv')

ggsave('female_plac_FINAL_volcano_plot_for_Alexis_paper.svg')





#### reading in female brain

m_brain <- read.csv('pearson_corr_w_triglycerides_f_brain.csv')

m_brain$log10_pval <- -log10(m_brain$pval)


### binning the correlation coefficients into 100 bins
bins <- bin_data(m_brain$CorrCoef, bins = 50, returnDT = TRUE)

### creating a new dataframe out of these bins which assigns a numeric value to each bin
binz <- unique(bins$Bin)
final_bins <- data.frame(binz)
final_bins$index <- seq(1, 50, by=1)


bins_final <- merge(bins, final_bins, by.x = 'Bin', by.y = 'binz')
### merging the bin dataframe with the final dataframe
new_df <- merge(m_brain, bins_final, by.x = 'CorrCoef', by.y = 'BinVal')
##drop the duplicates from this dataframe bc there is a weird error i cannot figure out
new_df <- new_df[!duplicated(new_df$Gene), ]


avg <- new_df %>%
  group_by(index) %>%
  summarize(binned_corr_coef = mean(CorrCoef))

final_merged <- merge(new_df, avg, by = 'index')

ggplot(data = final_merged, aes(x = binned_corr_coef, y =log10_pval, label = Gene))+
  geom_quasirandom(data = subset(final_merged, log10_pval > -log10(.05)), size = 2, color = 'red', width = .2, dodge.width = .5)+
  geom_quasirandom(data = subset(final_merged, log10_pval < -log10(.05)), size = 1, color = 'black', width = .2, dodge.width = .5)+
  geom_hline(yintercept = -log10(.05), linetype = 'dotted', size = 1.5)+
  geom_label_repel(data = subset(final_merged, Gene == 'DCX' | Gene == 'SLC4A7' | Gene == 'KCNA6' | Gene == 'KCNH7' | 
                                   Gene == 'TGFB3' | Gene == 'TGFBR2' | Gene == 'SMAD3' | 
                                   Gene == 'BMP4' | Gene == 'CD44' 
  ), box.padding = 2, max.overlaps = Inf)+
  #geom_label_repel(data = subset(final_merged, log10_pval > 3 & CorrCoef > .4 | CorrCoef < -.5),
  #box.padding   = 0.35, 
  #point.padding = 0.5,
  #segment.color = 'grey50')+
  ggtitle('Female Brain')+
  xlab('Binned CorrCoef')+
  ylab('Log10_Pval')+
  theme_bw()+
  scale_y_continuous(limits = c(0,6), expand = c(0, 0))+
  theme(plot.title = element_text(hjust = 0.5))

write.csv(final_merged, 'female_brain_w_bins.csv')

ggsave('female_brain_FINAL_volcano_plot_for_Alexis_paper.svg')

