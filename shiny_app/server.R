

######################################
###                                ###
###     server.R file for          ###
###       human-fetal-seq          ###
###        application             ###
###                                ###
######################################



server <- function (input, output, session) {
  
  updateSelectizeInput(session, 'entered_genes', choices = choice_genes, server = TRUE)
  
  output$landing_text <- renderUI ({
    
    x <- list(tags$h2("Welcome to Human Fetal RNA Seq!", style = 'text-align: center'),
              tags$h5("This website contains interactive visualizations and statistics of our bulk RNA sequencing data from
                      matched human fetal brain and placenta tissue (72-82 days post conception).",
                      style='text-align:center'),
              tags$h3(
                tags$a(href = 'http://microglia-seq.vm.duke.edu/microglia-seq/human-fetal-RNASeq/fetal_seq_mobile/', 'MOBILE SITE HERE'), style='text-align: center'),
              img(src="microglia.gif", height = '400px', width='505px', style = "display: block; margin-left: auto; margin-right: auto;"),
              tags$h5(""),
              fluidRow(column(width = 6,tags$br(), img(src="FetalTissue.png", height='100%', width='100%')),
              column(width = 6, tags$br(), tags$h5("Under the", tags$strong("Plots"), "menu on the sidebar:"), tags$h5("The", tags$strong("Correlation Plots"), "tab will show you if your gene of interest is correlated with maternal
              triglyceride accumulation (a proxy for dietary fat intake)."),
              tags$h5(""),
              tags$h5("The", tags$strong("Sex Differences"), "tab will show you if your gene of interest is differentially
                      expressed by sex in either brain or placenta."),
              tags$h5(""),
              tags$h5("The", tags$strong("Data Table"), "tab will provide the TPM values for your gene of interest alongside available data (sex,
              5-HT levels, maternal triglyceride accumulation) for each tissue."))),
              tags$h5(),
              tags$h3("",
                      "Check out the preprint",
                      tags$a(href="https://www.biorxiv.org/content/10.1101/2021.11.12.468408v2", 
                             "here!"), style = 'text-align: center')
              )
    tagList(x)
  })
  
  output$background_text <- renderUI({
    aa <- list(tags$h1("The Backstory", style='text-align: center'),
               tags$h4(),
               tags$h3("The",
                       tags$a(href="http://bilbolab.com/", target = '_blank', "Bilbo Lab"),
                       "is dedicated to understanding sex differences in developmental neuroimmune interactions", style = 'text-align: center'),
               tags$h3(""),
               tags$h4(
      HTML(paste("Recently, we uncovered that maternal high fat diet leads to endotoxin accumulation
              in fetal tissue in mice, resulting in increased Tlr4 signaling in fetal tissues.", 
                 tags$span(style="color:#05575D", " In male offspring"),"this leads to increased Tlr4-dependent
                           microglial phagocytosis of serotonin (5-HT) neurons in the fetal dorsal raphe nucleus (DRN), 
                           persistently decreased 5-HT in adult offspring brains, and diminished reward-seeking behavior.
                           Microglia in", tags$span(style="color:purple", "female offspring"), "show signs of increased
                           phagocytic activity, but the target CNS cell population is unknown.")), style='text-align:justify'),
      
      
      ##alexis added this 7/22 to try and adjust aspect ratios on mobile! needs testing. had to add a bunch of breaks after because otherwise it 
      ## was text wrapping which is weird. also can't get it centered - tried lots of things (margin: auto,etc)
      
      div(class ="col-sm-8", 
      img(src="humanVSmouse.png", height = '100px', width='460px', style = "display: block; margin-left: auto; margin-right: auto;")),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$h4(
        HTML(paste("We were interested in determining to what degree these findings could be translated to human populations, 
              and we obtained matched human fetal brain, placenta, and decidua (maternal placenta) from elective terminations
              matching the developmental timeline we followed in mice (embryonic days 14-16 in mice is roughly equivalent
              to 63-81 days post conception in humans). Given that no clinical data were available from these tissues, we assessed
              decidual triglyceride accumulation as a proxy for maternal dietary fat intake, and validated this method in our mouse model.
              We then correlated fetal tissue 5-HT levels with maternal triglyceride accumulation and found that fetal brain 5-HT was
              significantly negatively correlated with maternal triglyceride accumulation", tags$span(style="color:#05575D", "in males only."),
                   "We then performed bulk RNA-
              sequencing on the matched fetal brain and placenta tissues.")), style='text-align:justify'),
      img(src="summary.png", height="290px", width="620px",style = "display: block; margin-left: auto; margin-right: auto;"))
    tagList(aa)
  })
  
  
  output$table1 <- renderDT(options = list(pageLength = 35),{
    
    gene_of_interest <- input$entered_genes
    
    goi_df_brain <- brain_tpm_df[brain_tpm_df$hgnc_symbol == gene_of_interest, ]
    goi_df_placenta <- placenta_tpm_df[placenta_tpm_df$hgnc_symbol == gene_of_interest, ]
    goi_df_brain_t <- data.frame(t(goi_df_brain))
    goi_df_placenta_t <- data.frame(t(goi_df_placenta))
    goi_df_brain_t$SampleID <- rownames(goi_df_brain_t)
    goi_df_placenta_t$SampleID <- rownames(goi_df_placenta_t)
    for (row in rownames(goi_df_brain_t)) {
      goi_df_brain_t$SampleID[goi_df_brain_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
    }
    for (row in rownames(goi_df_placenta_t)) { 
      goi_df_placenta_t$SampleID[goi_df_placenta_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
    }
    brain_merged <- merge(goi_df_brain_t, metadata, by.x="SampleID", by.y="sampleID")
    brain_merged <- brain_merged %>%
      rename(Brain = starts_with("X"),
             Sex = sex,
             Age = age,
             Maternal_Triglycerides = trig,
             Brain_5HT = brain_5ht,
             Placenta_5HT = plac_5ht)
    
    placenta_merged <- merge(goi_df_placenta_t, metadata, by.x="SampleID", by.y="sampleID")
    placenta_merged <- placenta_merged %>%
      rename(Placenta = starts_with("X"),
             Sex = sex,
             Age = age,
             Maternal_Triglycerides = trig,
             Brain_5HT = brain_5ht,
             Placenta_5HT = plac_5ht)
    
    
    both_merged <- merge(placenta_merged, goi_df_brain_t, by.x = "SampleID", by.y = "SampleID", all.x = TRUE, all.y = TRUE)
    both_merged <- both_merged %>%
      rename(Brain = starts_with("X"))
    
    merged <- both_merged[order(both_merged$Sex),]
    merged_final <- subset(merged, select=-c(exclude, Notes))
   merged_final
  })
  
  
  
##Ben's original 
##  output$plot1 <- renderPlotly({
##    
##    goi <- input$entered_genes
##    goi_df <- brain_df[brain_df$hgnc_symbol == goi,]
##    goi_df <- data.frame(t(goi_df))
##    goi_df$sampleIDs <- rownames(goi_df)
##    rownames(goi_df) <- NULL
##    goi_df <- goi_df[6:length(rownames(goi_df)) - 1,]
##    colnames(goi_df)[1] <- goi
##    goi_df$sampleID <- goi_df$sampleIDs
##    goi_df$region <- goi_df$sampleIDs
##    for (sample in goi_df$sampleIDs) {
##      goi_df$sampleID[goi_df$sampleIDs == sample] <- str_split(sample, '_')[[1]][1]
##      goi_df$region[goi_df$sampleIDs == sample] <- str_split(sample, '_')[[1]][2]
##    }
##  
##  merged <- merge(goi_df, metadata, by.x="sampleID", by.y="sampleID")
##  merged[, goi] <- sapply(merged[, goi], as.numeric)
##  merged$sex <- factor(merged$sex, levels = c('Male', 'Female'))

## p <- ggscatter(merged, x = 'trig', y = goi, add = 'reg.line', conf.int = TRUE, cor.method = "pearson")+
##      ggtitle(paste(goi, 'Brain', sep = ' '))+
##      facet_wrap(~sex)+
##      corr_plot_settings
      
  
##  fig <- ggplotly(p, height = shinybrowser::get_height()-500, width = shinybrowser::get_width()-500)
##  fig %>% layout(margin = list(b = 90))
  
##  })
  
  output$plot1 <- renderPlotly({
    gene_of_interest <- input$entered_genes
    
    goi_df_brain <- brain_tpm_df[brain_tpm_df$hgnc_symbol == gene_of_interest, ]
    goi_df_placenta <- placenta_tpm_df[placenta_tpm_df$hgnc_symbol == gene_of_interest, ]
    
    if (dim(goi_df_placenta)[1]==0 && dim(goi_df_brain)[1]==0) {
      print("This gene is not present in the dataset")
      rm(goi_df_placenta)
      rm(goi_df_brain)
    }else if (dim(goi_df_placenta)[1]==0) {
      print("This gene is not present in the placenta")
      rm(goi_df_placenta)
      goi_df_brain_t <- data.frame(t(goi_df_brain))
      goi_df_brain_t$SampleID <- rownames(goi_df_brain_t)
      for (row in rownames(goi_df_brain_t)) {
        goi_df_brain_t$SampleID[goi_df_brain_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
        
      }
      brain_merged <- merge(goi_df_brain_t, metadata, by.x="SampleID", by.y="sampleID")
      brain_merged <- brain_merged %>%
        rename(Brain = starts_with("X"),
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      brain_merged$Brain <- as.numeric(brain_merged$Brain)
      male_brain_data <- subset(brain_merged, Sex == "Male")
      female_brain_data <- subset(brain_merged, Sex == "Female")
      malebrainStat <- cor.test(male_brain_data$Maternal_Triglycerides, male_brain_data$Brain)
      malebrainpval <- round(malebrainStat$p.value,3)
      malebrainpval <- paste("p=",malebrainpval)
      malebrainR <- round(malebrainStat$estimate, 2)
      malebrainR <- as.numeric(malebrainR)
      malebrainR <- paste("R=", malebrainR)
      malegraphstat <- as.character(paste(malebrainR, malebrainpval))
      malebrain_TPM_max = max(male_brain_data$Brain)
      femalebrainStat <- cor.test(female_brain_data$Maternal_Triglycerides, female_brain_data$Brain)
      femalebrainpval <- round(femalebrainStat$p.value,3)
      femalebrainpval <- paste("p=",femalebrainpval)
      femalebrainR <- round(femalebrainStat$estimate, 2)
      femalebrainR <- as.numeric(femalebrainR)
      femalebrainR <- paste("R=", femalebrainR)
      femalegraphstat <- as.character(paste(femalebrainR, femalebrainpval))
      femalebrain_TPM_max = max(female_brain_data$Brain)
      Male_Brain_Corr <- ggplotly(ggplot(male_brain_data, aes(x = Maternal_Triglycerides, y = Brain))+
                                    geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                                fill = "lightcyan3", alpha=0.5)+
                                    geom_point(color="turquoise4", alpha=0.5)+
                                    geom_line(stat="smooth", method="lm",color = "lightcyan4",
                                              size=.5)+
                                    ggtitle(label = gene_of_interest)+
                                    geom_text(x=.8, y=malebrain_TPM_max, aes(label=paste(malebrainR,",",malebrainpval)), size=3)+
                                    corr_plot_settings)
      
      Female_Brain_Corr <- ggplotly(ggplot(female_brain_data, aes(x = Maternal_Triglycerides, y = Brain))+
                                      geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                                  fill = "mistyrose3", alpha=0.5)+
                                      geom_point(color="rosybrown3", alpha=0.5)+
                                      geom_line(stat="smooth", method="lm",color = "rosybrown",
                                                size=.5)+
                                      ggtitle(label = gene_of_interest)+
                                      geom_text(x=1.1, y=femalebrain_TPM_max, aes(label=paste(femalebrainR,",",femalebrainpval)), size=3)+
                                      corr_plot_settings)
      splitCorrLayout <- subplot(Female_Brain_Corr, Male_Brain_Corr, titleY=TRUE, titleX = TRUE,
                                 margin=0.1)
      
      annotations = list( 
        list( 
          x = 0.2,  
          y = 1.0, 
          text = "Female Brain",
          xref = "paper",  
          yref = "paper",  
          xanchor = "center",  
          yanchor = "bottom",  
          showarrow = FALSE 
        ),  
        list( 
          x = 0.8,  
          y = 1,  
          text = "Male Brain",  
          xref = "paper",  
          yref = "paper",  
          xanchor = "center",  
          yanchor = "bottom",  
          showarrow = FALSE 
        ))
      ggplotly(splitCorrLayout) %>% layout(annotations = annotations)
    } else if (dim(goi_df_brain)[1]==0) {
      print("This gene is not present in the brain")
      rm(goi_df_brain)
      goi_df_placenta_t <- data.frame(t(goi_df_placenta))
      goi_df_placenta_t$SampleID <- rownames(goi_df_placenta_t)
      for (row in rownames(goi_df_placenta_t)) { 
        goi_df_placenta_t$SampleID[goi_df_placenta_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
        
      }
      placenta_merged <- merge(goi_df_placenta_t, metadata, by.x="SampleID", by.y="sampleID")
      placenta_merged <- placenta_merged %>%
        rename(Placenta = starts_with("X"),
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      placenta_merged$Placenta <- as.numeric(placenta_merged$Placenta)
      male_plac_data <- subset(placenta_merged, Sex == "Male")
      female_plac_data <- subset(placenta_merged, Sex == "Female")
      malePlacentaStat <- cor.test(male_plac_data$Maternal_Triglycerides, male_plac_data$Placenta)
      malePlacentapval <- round(malePlacentaStat$p.value,3)
      malePlacentapval <- paste("p=",malePlacentapval)
      malePlacentaR <- round(malePlacentaStat$estimate, 2)
      malePlacentaR <- as.numeric(malePlacentaR)
      malePlacentaR <- paste("R=", malePlacentaR)
      malegraphstat <- as.character(paste(malePlacentaR, malePlacentapval))
      malePlacenta_TPM_max = max(male_plac_data$Placenta)
      femalePlacentaStat <- cor.test(female_plac_data$Maternal_Triglycerides, female_plac_data$Placenta)
      femalePlacentapval <- round(femalePlacentaStat$p.value,3)
      femalePlacentapval <- paste("p=",femalePlacentapval)
      femalePlacentaR <- round(femalePlacentaStat$estimate, 2)
      femalePlacentaR <- as.numeric(femalePlacentaR)
      femalePlacentaR <- paste("R=", femalePlacentaR)
      femalegraphstat <- as.character(paste(femalePlacentaR, femalePlacentapval))
      femalePlacenta_TPM_max = max(female_plac_data$Placenta)
      Male_Plac_Corr <- ggplotly(ggplot(male_plac_data, aes(x = Maternal_Triglycerides, y = Placenta))+
                                   geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                               fill = "lightcyan3", alpha=0.5)+
                                   geom_point(color="turquoise4", alpha=0.5)+
                                   geom_line(stat="smooth", method="lm",color = "lightcyan4",
                                             size=.5)+
                                   ggtitle(label = gene_of_interest)+
                                   geom_text(x=.8, y=malePlacenta_TPM_max, aes(label=paste(malePlacentaR,",",malePlacentapval)), size=3)+
                                   corr_plot_settings)
      
      Female_Plac_Corr <- ggplotly(ggplot(female_plac_data, aes(x = Maternal_Triglycerides, y = Placenta))+
                                     geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                                 fill = "mistyrose3", alpha=0.5)+
                                     geom_point(color="rosybrown3", alpha=0.5)+
                                     geom_line(stat="smooth", method="lm",color = "rosybrown",
                                               size=.5)+
                                     ggtitle(label = gene_of_interest)+
                                     geom_text(x=1.1, y=femalePlacenta_TPM_max, aes(label=paste(femalePlacentaR,",",femalePlacentapval)), size=3)+
                                     corr_plot_settings)
      
      splitCorrLayout <- subplot(Female_Plac_Corr, Male_Plac_Corr, titleY=TRUE, titleX = TRUE,
                                 margin = .1) 
      
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
        ))
      ggplotly(splitCorrLayout) %>% layout(annotations = annotations)
    } else {
      
      
      
      goi_df_brain_t <- data.frame(t(goi_df_brain))
      goi_df_placenta_t <- data.frame(t(goi_df_placenta))
      
      
      ## rownames(goi_df_brain_t)
      
      goi_df_brain_t$SampleID <- rownames(goi_df_brain_t)
      goi_df_placenta_t$SampleID <- rownames(goi_df_placenta_t)
      
      
      
      for (row in rownames(goi_df_brain_t)) {
        goi_df_brain_t$SampleID[goi_df_brain_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
        
      }
      
      
      for (row in rownames(goi_df_placenta_t)) { 
        goi_df_placenta_t$SampleID[goi_df_placenta_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
        
      }
      
      
      
      brain_merged <- merge(goi_df_brain_t, metadata, by.x="SampleID", by.y="sampleID")
      brain_merged <- brain_merged %>%
        rename(Brain = starts_with("X"),
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      
      placenta_merged <- merge(goi_df_placenta_t, metadata, by.x="SampleID", by.y="sampleID")
      placenta_merged <- placenta_merged %>%
        rename(Placenta = starts_with("X"),
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      
      brain_merged$Brain <- as.numeric(brain_merged$Brain)
      placenta_merged$Placenta <- as.numeric(placenta_merged$Placenta)
      
      male_brain_data <- subset(brain_merged, Sex == "Male")
      male_plac_data <- subset(placenta_merged, Sex == "Male")
      female_brain_data <- subset(brain_merged, Sex == "Female")
      female_plac_data <- subset(placenta_merged, Sex == "Female")
      
      malebrainStat <- cor.test(male_brain_data$Maternal_Triglycerides, male_brain_data$Brain)
      malebrainpval <- round(malebrainStat$p.value,3)
      malebrainpval <- paste("p=",malebrainpval)
      malebrainR <- round(malebrainStat$estimate, 2)
      malebrainR <- as.numeric(malebrainR)
      malebrainR <- paste("R=", malebrainR)
      malegraphstat <- as.character(paste(malebrainR, malebrainpval))
      malebrain_TPM_max = max(male_brain_data$Brain)
      
      femalebrainStat <- cor.test(female_brain_data$Maternal_Triglycerides, female_brain_data$Brain)
      femalebrainpval <- round(femalebrainStat$p.value,3)
      femalebrainpval <- paste("p=",femalebrainpval)
      femalebrainR <- round(femalebrainStat$estimate, 2)
      femalebrainR <- as.numeric(femalebrainR)
      femalebrainR <- paste("R=", femalebrainR)
      femalegraphstat <- as.character(paste(femalebrainR, femalebrainpval))
      femalebrain_TPM_max = max(female_brain_data$Brain)
      
      malePlacentaStat <- cor.test(male_plac_data$Maternal_Triglycerides, male_plac_data$Placenta)
      malePlacentapval <- round(malePlacentaStat$p.value,3)
      malePlacentapval <- paste("p=",malePlacentapval)
      malePlacentaR <- round(malePlacentaStat$estimate, 2)
      malePlacentaR <- as.numeric(malePlacentaR)
      malePlacentaR <- paste("R=", malePlacentaR)
      malegraphstat <- as.character(paste(malePlacentaR, malePlacentapval))
      malePlacenta_TPM_max = max(male_plac_data$Placenta)
      
      femalePlacentaStat <- cor.test(female_plac_data$Maternal_Triglycerides, female_plac_data$Placenta)
      femalePlacentapval <- round(femalePlacentaStat$p.value,3)
      femalePlacentapval <- paste("p=",femalePlacentapval)
      femalePlacentaR <- round(femalePlacentaStat$estimate, 2)
      femalePlacentaR <- as.numeric(femalePlacentaR)
      femalePlacentaR <- paste("R=", femalePlacentaR)
      femalegraphstat <- as.character(paste(femalePlacentaR, femalePlacentapval))
      femalePlacenta_TPM_max = max(female_plac_data$Placenta)
      
      
      Male_Brain_Corr <- ggplotly(ggplot(male_brain_data, aes(x = Maternal_Triglycerides, y = Brain))+
                                    geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                                fill = "lightcyan3", alpha=0.5)+
                                    geom_point(color="turquoise4", alpha=0.5)+
                                    geom_line(stat="smooth", method="lm",color = "lightcyan4",
                                              size=.5)+
                                    ggtitle(label = gene_of_interest)+
                                    geom_text(x=.8, y=malebrain_TPM_max, aes(label=paste(malebrainR,",",malebrainpval)), size=3)+
                                    corr_plot_settings)
      
      Female_Brain_Corr <- ggplotly(ggplot(female_brain_data, aes(x = Maternal_Triglycerides, y = Brain))+
                                      geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                                  fill = "mistyrose3", alpha=0.5)+
                                      geom_point(color="rosybrown3", alpha=0.5)+
                                      geom_line(stat="smooth", method="lm",color = "rosybrown",
                                                size=.5)+
                                      geom_text(x=1.1, y=femalebrain_TPM_max, aes(label=paste(femalebrainR,",",femalebrainpval)), size=3)+
                                      corr_plot_settings)
      
      Male_Plac_Corr <- ggplotly(ggplot(male_plac_data, aes(x = Maternal_Triglycerides, y = Placenta))+
                                   geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                               fill = "lightcyan3", alpha=0.5)+
                                   geom_point(color="turquoise4", alpha=0.5)+
                                   geom_line(stat="smooth", method="lm",color = "lightcyan4",
                                             size=.5)+
                                   geom_text(x=.8, y=malePlacenta_TPM_max, aes(label=paste(malePlacentaR,",",malePlacentapval)), size=3)+
                                   corr_plot_settings)
      
      Female_Plac_Corr <- ggplotly(ggplot(female_plac_data, aes(x = Maternal_Triglycerides, y = Placenta))+
                                     geom_ribbon(stat="smooth", method = "lm", se=TRUE,
                                                 fill = "mistyrose3", alpha=0.5)+
                                     geom_point(color="rosybrown3", alpha=0.5)+
                                     geom_line(stat="smooth", method="lm",color = "rosybrown",
                                               size=.5)+
                                     geom_text(x=1.1, y=femalePlacenta_TPM_max, aes(label=paste(femalePlacentaR,",",femalePlacentapval)), size=3)+
                                     corr_plot_settings)
      brplot <- ggplotly(ggplot())
      
      
      
      splitCorrLayout <- subplot(Female_Plac_Corr, brplot, Male_Plac_Corr, brplot, brplot, brplot, Female_Brain_Corr, brplot, Male_Brain_Corr, nrows=3, 
                                 titleY=TRUE, titleX = TRUE, margin=0.04, heights = c(0.45, 0.1, 0.45), widths = c(0.45, 0.1, 0.45))
      
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
          x = 0.75,  
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
          x = 0.75,  
          y = 0.4,
          text = "Male Brain",  
          xref = "paper",  
          yref = "paper",  
          xanchor = "center",  
          yanchor = "bottom",
          showarrow = FALSE 
        ))
      
      
      ggplotly(splitCorrLayout) %>% layout(annotations = annotations, autosize = F, width = 800, height = 750)
    }
  })
    
## Next is sex differences plot  
  
  output$sex_diff_plot <- renderPlotly({
    
    
    gene_of_interest <- input$entered_genes
    
    
    goi_df_brain <- brain_tpm_df[brain_tpm_df$hgnc_symbol == gene_of_interest, ]
    goi_df_placenta <- placenta_tpm_df[placenta_tpm_df$hgnc_symbol == gene_of_interest, ]
    
    if (dim(goi_df_placenta)[1]==0 && dim(goi_df_brain)[1]==0) {
      print("This gene is not present in the dataset")
      rm(goi_df_placenta)
      rm(goi_df_brain)
    }else if (dim(goi_df_placenta)[1]==0) {
      print("This gene is not present in the placenta")
      rm(goi_df_placenta)
      goi_df_brain_t <- data.frame(t(goi_df_brain))
      goi_df_brain_t$SampleID <- rownames(goi_df_brain_t)
      for (row in rownames(goi_df_brain_t)) {
        goi_df_brain_t$SampleID[goi_df_brain_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
      }
      brain_merged <- merge(goi_df_brain_t, metadata, by.x="SampleID", by.y="sampleID")
      brain_merged <- brain_merged %>%
        rename(Brain = starts_with("X"), 
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      brain_merged$Brain <- as.numeric(brain_merged$Brain)
      brain_merged_grouped <- brain_merged %>%
        group_by(Sex) %>%
        summarise(sem_brain = sem(Brain), Brain = mean(Brain))
      brain_pval <- 
        t.test(brain_merged$Brain[brain_merged$Sex == "Male"],
               brain_merged$Brain[brain_merged$Sex == "Female"])
      brain_pval_graph <- round(brain_pval$p.value, 4)
      brain_pval_graph <- paste("p=",brain_pval_graph)
      b_max <- max(brain_merged$Brain)
      b_min <- min(brain_merged$Brain)
      brainBox <- ggplotly(ggplot(brain_merged, aes(x = Sex, y = Brain, fill = Sex))+
                             geom_boxplot(outlier.color = NA, outlier.size = 0, outlier.shape = NA)+
                             coord_cartesian(ylim = c(b_min, b_max))+
                             scale_fill_manual(values = c("mistyrose3", "lightcyan3"))+
                             geom_jitter(data = brain_merged, width = 0.1, alpha = 0.5)+
                             ggtitle(paste(label = gene_of_interest,"<br>Brain"))+
                             geom_text(aes(x=1.5, y=b_max, label=as.character(brain_pval_graph)), size=3)+
                             box_plot_theme)
      brainBox$x$data[[1]]$marker$opacity = 0 
      brainBox$x$data[[2]]$marker$opacity = 0
      ggplotly(brainBox) %>% layout(yaxis = yy, width = 500)
    } else if (dim(goi_df_brain)[1]==0) {
      print("This gene is not present in the brain")
      rm(goi_df_brain)
      goi_df_placenta_t <- data.frame(t(goi_df_placenta))
      goi_df_placenta_t$SampleID <- rownames(goi_df_placenta_t)
      for (row in rownames(goi_df_placenta_t)) {
        goi_df_placenta_t$SampleID[goi_df_placenta_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
        
      }
      placenta_merged <- merge(goi_df_placenta_t, metadata, by.x="SampleID", by.y="sampleID")
      placenta_merged <- placenta_merged %>%
        rename(Placenta = starts_with("X"), 
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      placenta_merged$Placenta <- as.numeric(placenta_merged$Placenta)
      placenta_merged_grouped <- placenta_merged %>%
        group_by(Sex) %>%
        summarise(sem_placenta = sem(Placenta), Placenta = mean(Placenta))
      placenta_pval <- 
        t.test(placenta_merged$Placenta[placenta_merged$Sex == "Male"],
               placenta_merged$Placenta[placenta_merged$Sex == "Female"])
      placenta_pval_graph <- round(placenta_pval$p.value, 4)
      placenta_pval_graph <- paste("p=",placenta_pval_graph)
      p_max <- max(placenta_merged$Placenta)
      p_min <- min(placenta_merged$Placenta)
      placentaBox <- ggplotly(ggplot(placenta_merged, aes(x = Sex, y = Placenta, fill = Sex))+
                                geom_boxplot(outlier.shape = NA, outlier.size = NA, outlier.color = NA)+
                                coord_cartesian(ylim = c(p_min, p_max))+
                                scale_fill_manual(values = c("mistyrose3", "lightcyan3"))+
                                geom_jitter(data = placenta_merged, width = 0.1, alpha = 0.5)+
                                ggtitle(paste(label = gene_of_interest,"<br>Placenta"))+
                                geom_text(aes(x=1.5, y=p_max, label=as.character(placenta_pval_graph)), size=3)+
                                box_plot_theme)
      placentaBox$x$data[[1]]$marker$opacity = 0 
      placentaBox$x$data[[2]]$marker$opacity = 0
      ggplotly(placentaBox) %>% layout(yaxis = yy, width = 500)
    } else {
      
      goi_df_brain_t <- data.frame(t(goi_df_brain))
      goi_df_placenta_t <- data.frame(t(goi_df_placenta))
      
      ## this returns a list not a dataframe. we needed it to be a dataframe so we tell it to be a dataframe
      ## then we need to make the row names into an actual row b/c they behave different.
      ## to check, do 
      
      ## rownames(goi_df_brain_t)
      
      goi_df_brain_t$SampleID <- rownames(goi_df_brain_t)
      goi_df_placenta_t$SampleID <- rownames(goi_df_placenta_t)
      
      
      
      for (row in rownames(goi_df_brain_t)) {
        goi_df_brain_t$SampleID[goi_df_brain_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
      }
      
      
      for (row in rownames(goi_df_placenta_t)) { 
        goi_df_placenta_t$SampleID[goi_df_placenta_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
      }
      
      
      ## create a new dataframe w/the metadata plus a column by merging metadata w/TPM data 
      
      brain_merged <- merge(goi_df_brain_t, metadata, by.x="SampleID", by.y="sampleID")
      brain_merged <- brain_merged %>%
        rename(Brain = starts_with("X"))
      
      both_merged <- merge(goi_df_placenta_t, brain_merged, by.x = "SampleID", by.y = "SampleID")
      both_merged <- both_merged %>%
        rename(Placenta = starts_with("X"),
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      
      brain_merged <- brain_merged %>%
        rename(Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      
      placenta_merged <- merge(goi_df_placenta_t, metadata, by.x="SampleID", by.y="sampleID")
      placenta_merged <- placenta_merged %>%
        rename(Placenta = starts_with("X"),
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      
      brain_merged$Brain <- as.numeric(brain_merged$Brain)
      placenta_merged$Placenta <- as.numeric(placenta_merged$Placenta)
      
      brain_merged_grouped <- brain_merged %>%
        group_by(Sex) %>%
        summarise(sem_brain = sem(Brain), Brain = mean(Brain))
      
      placenta_merged_grouped <- placenta_merged %>%
        group_by(Sex) %>%
        summarise(sem_placenta = sem(Placenta), Placenta = mean(Placenta))
      
      brain_pval <- 
        t.test(brain_merged$Brain[brain_merged$Sex == "Male"],
               brain_merged$Brain[brain_merged$Sex == "Female"])
      brain_pval_graph <- round(brain_pval$p.value, 4)
      brain_pval_graph <- paste("p=",brain_pval_graph)
      
      placenta_pval <- 
        t.test(placenta_merged$Placenta[placenta_merged$Sex == "Male"],
               placenta_merged$Placenta[placenta_merged$Sex == "Female"])
      placenta_pval_graph <- round(placenta_pval$p.value,4)
      placenta_pval_graph <- paste("p=",placenta_pval_graph)
      
      
      b_max <- max(brain_merged$Brain)
      p_max <- max(placenta_merged$Placenta)
      combined <- data.frame(b_max, p_max)
      y_max <- max(combined)
      
      
      b_min <- min(brain_merged$Brain)
      p_min <- min(placenta_merged$Placenta)
      combined_min <- data.frame(b_min, p_min)
      y_min <- min(combined_min)
      
      
      brainBox <- ggplotly(ggplot(brain_merged, aes(x = Sex, y = Brain, fill = Sex))+
                             geom_boxplot(outlier.color = NA, outlier.size = 0, outlier.shape = NA)+
                             coord_cartesian(ylim = c(y_min, y_max))+
                             scale_fill_manual(values = c("mistyrose3", "lightcyan3"))+
                             geom_jitter(data = brain_merged, width = 0.1, alpha = 0.5)+
                             ggtitle(label = gene_of_interest)+
                             geom_text(aes(x=1.5, y=y_max, label=as.character(brain_pval_graph)), size=3)+
                             box_plot_theme)
      
      placentaBox <- ggplotly(ggplot(placenta_merged, aes(x = Sex, y = Placenta, fill = Sex))+
                                geom_boxplot(outlier.color = NA, outlier.size = 0, outlier.shape = NA)+
                                coord_cartesian(ylim = c(y_min, y_max))+
                                scale_fill_manual(values = c("mistyrose3", "lightcyan3"))+
                                geom_jitter(data = placenta_merged, width = 0.1, alpha = 0.5)+
                                geom_text(aes(x=1.5, y=y_max, label=as.character(placenta_pval_graph)), size=3)+
                                box_plot_theme)
      
      ## setting outlier points to zero
      
      brainBox$x$data[[1]]$marker$opacity = 0 
      placentaBox$x$data[[1]]$marker$opacity = 0 
      brainBox$x$data[[2]]$marker$opacity = 0 
      placentaBox$x$data[[2]]$marker$opacity = 0 
      
      brplot2 <- ggplotly(ggplot())
      
      
      boxLayout <- subplot(style(placentaBox, showlegend = FALSE), brplot2, brainBox, nrows=1, widths = c(0.45, 0.1, 0.45),
                           margin = 0.04)
      
      annotations = list(
        list(
          x=0.2, 
          y=1, 
          text = "Placenta",
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "bottom",
          showarrow = FALSE
        ),
        list(
          x=0.8,
          y=1,
          text = "Brain",
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "bottom",
          showarrow = FALSE
        ))
      
      
      ggplotly(boxLayout) %>% layout(yaxis = yy, annotations = annotations)
      }
  })
  
  output$data_table_download <- downloadHandler(
    
    filename = function() {paste('data_output', '_', list(input$entered_genes), '.csv', sep = '')},
    content = function(file) {
      
      
      
      gene_of_interest <- input$entered_genes
      
      goi_df_brain <- brain_tpm_df[brain_tpm_df$hgnc_symbol == gene_of_interest, ]
      goi_df_placenta <- placenta_tpm_df[placenta_tpm_df$hgnc_symbol == gene_of_interest, ]
      goi_df_brain_t <- data.frame(t(goi_df_brain))
      goi_df_placenta_t <- data.frame(t(goi_df_placenta))
      goi_df_brain_t$SampleID <- rownames(goi_df_brain_t)
      goi_df_placenta_t$SampleID <- rownames(goi_df_placenta_t)
      for (row in rownames(goi_df_brain_t)) {
        goi_df_brain_t$SampleID[goi_df_brain_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
      }
      for (row in rownames(goi_df_placenta_t)) { 
        goi_df_placenta_t$SampleID[goi_df_placenta_t$SampleID == row] <- str_split(str_split(row, '_')[[1]][1], 'X')[[1]][2]
      }
      brain_merged <- merge(goi_df_brain_t, metadata, by.x="SampleID", by.y="sampleID")
      brain_merged <- brain_merged %>%
        rename(Brain = starts_with("X"),
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      
      placenta_merged <- merge(goi_df_placenta_t, metadata, by.x="SampleID", by.y="sampleID")
      placenta_merged <- placenta_merged %>%
        rename(Placenta = starts_with("X"),
               Sex = sex,
               Age = age,
               Maternal_Triglycerides = trig,
               Brain_5HT = brain_5ht,
               Placenta_5HT = plac_5ht)
      
      
      both_merged <- merge(goi_df_placenta_t, brain_merged, by.x = "SampleID", by.y = "SampleID")
      both_merged <- both_merged %>%
        rename(Placenta = starts_with("X"))
      
      merged <- both_merged[order(both_merged$Sex),]
      merged_final <- subset(merged, select=-c(exclude, Notes))
      merged_final
      
      write.csv(merged_final, file)
      
    })
  
  
  output$website_info_text <- renderUI({
    
    z <- list(
      fluidRow(
        column(width = 12,
               div(style = "height:50px;background-color: transparent;"))),
      tags$h3('All of the code used to analyze the data and produce this website can be found here: ', tags$a(
        href = 'https://github.com/bendevlin18/human-fetal-RNASeq', target = '_blank',
        tags$img(src="github-logo.png", title="GitHub", width="25", height="25"))),
      tags$h3("The raw and preprocessed data has been uploaded to the NCBI's GEO database", tags$a('GSE188872',
        href = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188872', target = '_blank')),
      tags$h4('If you run into any issues or have any questions about the website or the data contained within, contact us!'),
      fluidRow(
        column(width = 12,
               div(style = "height:50px;background-color: transparent;"))),
      fluidRow(
        column(width = 4, 
               tags$img(src="staci.jpg", title="Staci Image", width="225", height="231"),
               tags$p('Staci Bilbo, Ph.D'),
               tags$p('PI'),
               tags$a('staci.bilbo@duke.edu'),
               tags$a(
                 href = 'https://twitter.com/staci_bilbo', target = '_blank',
                 tags$img(src="twitter.svg", title="Twitter", width="25", height="25"))),
        column(width = 4, 
               tags$img(src="alexis.jpg", title="Staci Image", width="225", height="231"),
               tags$p('Alexis Ceasrine, Ph.D'),
               tags$p('Postdoc in the Bilbo Lab'),
               tags$a('alexis.ceasrine@duke.edu'),
               tags$a(
                 href = 'https://twitter.com/aceasrine', target = '_blank',
                 tags$img(src="twitter.svg", title="Twitter", width="25", height="25"))),
        column(width = 4, 
               tags$img(src="ben.png", title="Ben Image", width="225", height="225"),
               tags$p('Ben Devlin'),
               tags$p('Graduate Student in the Bilbo Lab'),
               tags$a('benjamin.devlin@duke.edu'),
               tags$a(
                 href = 'https://twitter.com/BenDevlin18', target = '_blank',
                 tags$img(src="twitter.svg", title="Twitter", width="25", height="25")))
               ),
      fluidRow(
        column(width = 12,
               div(style = "height:50px;background-color: transparent;")),
      tags$h3('Check out the Bilbo Lab Site!', tags$a('Lab Website',
        href = 'http://bilbolab.com', target = '_blank'))
        
      )
    )
    
    
    
    tagList(z)
    
    
    
  })
  
  
  
  
  
  
  
  
  
  
}