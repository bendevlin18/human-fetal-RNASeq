

rm(list = ls())

library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(plotly)
library(glue)

corr_plot_settings <- list(
  xlab('Log Maternal Triglycerides (mg/dL)'),
  ylab('TPM (Transcripts per million)'),
  theme_classic(),
  theme(rect = element_rect(fill = 'transparent'), 
        text = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11),
        plot.margin = margin(2, 2, 2, 2),
        axis.line = element_line(size = .75),
        axis.ticks = element_line(size = .6),
        axis.text.x = element_text(hjust = 1),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color = NA, col = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size=10),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 10),
        strip.background = element_rect(color="transparent", fill="transparent", linetype="solid")))

setwd('/Users/alexis/Library/CloudStorage/OneDrive-DukeUniversity/alexis_sequencing/shiny_app')


brain_tpm_df <- read.csv('brain_tpm_w_exclusions.csv')
placenta_tpm_df <- read.csv('placenta_tpm_w_exclusions.csv')



metadata <- read.csv('human_seq_metadata.csv')



gene_of_interest <- 'ADRB2'



goi_df_brain <- brain_tpm_df[brain_tpm_df$hgnc_symbol == gene_of_interest, ]
goi_df_placenta <- placenta_tpm_df[placenta_tpm_df$hgnc_symbol == gene_of_interest, ]



goi_df_brain_t <- data.frame(t(goi_df_brain))
goi_df_placenta_t <- data.frame(t(goi_df_placenta))


rownames(goi_df_brain_t)

goi_df_brain_t$SampleID <- rownames(goi_df_brain_t)
goi_df_placenta_t$SampleID <- rownames(goi_df_placenta_t)



for (row in rownames(goi_df_brain_t)) {
  goi_df_brain_t$SampleID[goi_df_brain_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
  
}


for (row in rownames(goi_df_placenta_t)) { 
  goi_df_placenta_t$SampleID[goi_df_placenta_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
  
}


for (sample in metadata$sampleID) {
  metadata$sampleID[metadata$sampleID == sample] <- str_split(sample, "-")[[1]][2]
}


brain_merged <- merge(goi_df_brain_t, metadata, by.x="SampleID", by.y="sampleID")
brain_merged <- brain_merged %>%
  rename(Brain = starts_with("X"))

placenta_merged <- merge(goi_df_placenta_t, metadata, by.x="SampleID", by.y="sampleID")
placenta_merged <- placenta_merged %>%
  rename(Placenta = starts_with("X"))


both_merged <- merge(goi_df_placenta_t, brain_merged, by.x = "SampleID", by.y = "SampleID")
both_merged <- both_merged %>%
  rename(Placenta = starts_with("X"),
         Sex = sex,
         Age = age,
         Maternal_Triglycerides = trig,
         Brain_5HT = brain_5ht,
         Placenta_5HT = plac_5ht)


both_merged$Brain <- as.numeric(both_merged$Brain)
both_merged$Placenta <- as.numeric(both_merged$Placenta)


##PlacentaPlots <- ggplotly(ggscatter(both_merged, x = 'Maternal_Triglycerides', y = 'Placenta', add = 'reg.line', conf.int = TRUE, cor.method = "pearson")+
##           facet_wrap(~Sex)+
##           xlab('Log Triglycerides (mg/dL)')+
##           ylab(paste(gene_of_interest, '(TPM)', sep = ' '))+
##           theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5))+
##           corr_plot_settings)

##BrainPlots <- ggplotly(ggscatter(both_merged, x = 'Maternal_Triglycerides', y = 'Brain', add = 'reg.line', conf.int = TRUE, cor.method = "pearson")+
##                            ggtitle(label = gene_of_interest)+
##                            facet_wrap(~Sex)+
##                            xlab('Log Triglycerides (mg/dL)')+
##                            ylab(paste(gene_of_interest, '(TPM)', sep = ' '))+
##                            theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5))+
##                            corr_plot_settings)

##CorrLayout <- subplot(PlacentaPlots, BrainPlots, nrows = 2)

##annotations = list(
##  list(
##    x=0.2, 
##    y=1, 
##    text = "Placenta",
    ##xref = "paper",
    ## yref = "paper",
##    xanchor = "center",
##    yanchor = "bottom",
##    showarrow = FALSE
##  ),
##  list(
##    x=0.8,
##    y=1,
##    text = "Brain",
##    xref = "paper",
##    yref = "paper",
##    xanchor = "center",
##    yanchor = "bottom",
##    showarrow = FALSE
##  ))



##ggplotly(CorrLayout) %>% layout(annotations = annotations)

male_data <- subset(both_merged, Sex == "Male")
female_data <- subset(both_merged, Sex == "Female")

Male_Brain_Corr <- ggplotly(ggplot(male_data, aes(x = Maternal_Triglycerides, y = Brain))+
                              geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                          fill = "lightcyan3", alpha=0.5)+
                              geom_point(color="turquoise4", alpha=0.5)+
                              geom_line(stat="smooth", method="lm",color = "lightcyan4",
                                        size=.5)+
##                                      cor.coef = TRUE,
##                                      cor.method = "pearson",
##                                      add.params = list(color = "blue", fill = "lightgray"))+
##)+
##                                      stat_cor(method = "pearson", label.x = 1, p.accuracy = 0.001, r.accuracy = 0.01)+
                                     ## cor.method = "pearson")+
                              ggtitle(label = gene_of_interest)+
                              corr_plot_settings+
  stat_cor(method = "pearson", label.x = 1, p.accuracy = 0.001, r.accuracy = 0.01))

ggplotly(Male_Brain_Corr)

Female_Brain_Corr <- ggplotly(ggscatter(female_data, x = 'Maternal_Triglycerides', y = 'Brain', 
                                        add = 'reg.line', conf.int = TRUE, cor.method = "pearson")+
                              corr_plot_settings)
Male_Plac_Corr <- ggplotly(ggscatter(male_data, x = 'Maternal_Triglycerides', y = 'Placenta', 
                                     add = 'reg.line', conf.int = TRUE, cor.method = "pearson")+
                              corr_plot_settings)
Female_Plac_Corr <- ggplotly(ggscatter(female_data, x = 'Maternal_Triglycerides', y = 'Placenta', 
                                       add = 'reg.line', conf.int = TRUE, cor.method = "pearson")+
                                corr_plot_settings)


splitCorrLayout <- subplot(Female_Plac_Corr, Male_Plac_Corr, Female_Brain_Corr, Male_Brain_Corr, nrows = 2, titleY=TRUE, titleX = TRUE,
                           margin=0.1)

annotations = list( 
  list( 
    x = 0.2,  
    y = 1.0,  
    text = "Female Placenta",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),  
  list( 
    x = 0.8,  
    y = 1,  
    text = "Male Placenta",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),  
  list( 
    x = 0.2,  
    y = 0.4,  
    text = "Female Brain",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),
  list( 
    x = 0.8,  
    y = 0.4,  
    text = "Male Brain",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ))



ggplotly(splitCorrLayout) %>% layout(annotations = annotations)




##plot(both_merged$Placenta, both_merged$Brain, main="Correlation of gene expression", xlab = "Placenta TPM", ylab = 
##       "Brain TPM")

##ggplot(both_merged, aes(x=Placenta, y=Brain, color = Maternal_Triglycerides))+
##  geom_point(size=3)


##ggscatter(both_merged, x = 'Placenta', y = 'Brain',
##          add = "reg.line", conf.int = TRUE,
##          cor.coef = TRUE, cor.method = "pearson" )
