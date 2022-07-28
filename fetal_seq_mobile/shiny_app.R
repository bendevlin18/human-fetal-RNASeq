

##### Load in necessary packages #####

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
brain_df <- read.csv('brain_tpm_w_exclusions.csv')
placenta_df <- read.csv('placenta_tpm_w_exclusions.csv')

brain_genes <- brain_df$hgnc_symbol
placenta_genes <- placenta_df$hgnc_symbol

z <- t(brain_df)
## batch through and clean colnames for the brain
for (cols in colnames(brain_df)) {
  
  if (grepl('tpm', cols, fixed=TRUE)) {
    
    colnames(brain_df)[colnames(brain_df) == cols] <- str_split(str_split(cols, 'X')[[1]][2], '_tpm')[[1]][1]
    
  }
}


z <- t(placenta_df)
## batch through and clean colnames for the placenta
for (cols in colnames(placenta_df)) {
  
  if (grepl('tpm', cols, fixed=TRUE)) {
    
    colnames(placenta_df)[colnames(placenta_df) == cols] <- str_split(str_split(cols, 'X')[[1]][2], '_tpm')[[1]][1]
    
  }
}



###### TISSUE AND SAMPLE METADATA #######

metadata <- read.csv('human_seq_metadata.csv')

### prepping metadata
metadata$sampleID <- metadata$ï..Sample.ID
for (sample in metadata$sampleID) {
  
  metadata$sampleID[metadata$sampleID == sample] <- str_split(sample, "-")[[1]][2]
  
}
metadata <- metadata[!grepl('y', metadata$exclude),]


#################################################




###### UI Section ######

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Hello!"),
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    ## input: select a gene
    textInput('entered_genes', label = 'Enter a gene', value = 'NEUROD1'),
    
    # Output: correlation plot
    plotOutput(outputId = "male_female_plot")
    
  )
)




####### Server section


server <- function(input, output, session) {
  
  goi = input$entered_genes
  
  goi_df <- brain_df[brain_df$hgnc_symbol %in% goi,]
  goi_df <- data.frame(t(goi_df))
  goi_df$sampleIDs <- rownames(goi_df)
  rownames(goi_df) <- NULL
  goi_df <- goi_df[4:81,]
  colnames(goi_df)[1] <- goi
  
  ## looping to extract the sampleID column and the region column
  goi_df$sampleID <- goi_df$sampleIDs
  goi_df$region <- goi_df$sampleIDs
  for (sample in goi_df$sampleIDs) {
    
    
    
    goi_df$sampleID[goi_df$sampleIDs == sample] <- str_split(str_split(sample, 'X')[[1]][2], '_')[[1]][1]
    
    goi_df$region[goi_df$sampleIDs == sample] <- str_split(str_split(sample, 'X')[[1]][2], '_')[[1]][2]
    
  }
  
  
  ## merging goi_df with the metadata
  merged <- merge(goi_df, metadata, by.x="sampleID", by.y="sampleID")
  merged <- transform(merged, `goi` = as.numeric(`goi`))
  merged$sex <- factor(merged$sex, levels = c('Male', 'Female'))
  
  output$male_female_plot <- renderPlot({
    ggscatter(merged, x = 'trig', y = 'TLR5', add = 'reg.line', conf.int = TRUE, cor.method = "pearson")+
      ggtitle(paste(gene_of_interest, 'Brain', sep = ' '))+
      facet_wrap(~sex)+
      xlab('Log Triglycerides (mg/dL)')+
      ylab(paste(gene_of_interest, '(TPM)', sep = ' '))+
      theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5))
    
  })
  
}







shinyApp(ui = ui, server = server)
