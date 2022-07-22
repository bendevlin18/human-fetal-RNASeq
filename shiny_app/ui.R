

######################################
###                                ###
###       ui.R file for            ###
###       human-fetal-seq          ###
###        application             ###
###                                ###
######################################

ui <-  dashboardPage(skin = 'blue',
                     dashboardHeader(title = 'Human Fetal RNASeq', titleWidth = 350),
                     dashboardSidebar(
                       sidebarMenu(
                         menuItem('Landing Page', tabName = 'landing', icon = icon('brain')),
                         menuItem('Background', tabName = 'background', icon = icon('file-text-o')),
                         menuItem('Plots', tabName = 'corr_plots', icon = icon('chart-line'), startExpanded = FALSE,
                                  selectizeInput('entered_genes', multiple = TRUE, choices = NULL, label = 'Search for gene of interest',
                                                 options = list(maxItems = 1, placeholder = 'Example: Cx3cr1')),
                         menuSubItem('Correlation Plots', tabName = 'br_corr_plots', icon = icon('chart-line')),
                         menuSubItem('Sex Differences', tabName = 'sex_diff', icon = icon('chart-bar')),
                         menuSubItem('Data Table', tabName = 'br_data_table', icon = icon('table'))
                         
                         )
                        
                                  )
                       
                     ),
                     dashboardBody(
                       shinybrowser::detect(),
                       
                       
            
                       tabItems(
                         tabItem('landing', uiOutput('landing_text')),
                         tabItem('background', uiOutput('background_text')),
                         tabItem('br_data_table', DTOutput("table1"),
                                 p(class = 'text-center', downloadButton('data_table_download', 'Download Data Table'))
                                 ),
                         tabItem('br_corr_plots',
                                 column(width = 12, align='center', withSpinner(plotlyOutput('plot1')))
                                 ),
                         tabItem('sex_diff', column(width = 12, align='center', withSpinner(plotlyOutput('sex_diff_plot')))
                         )
                       )))
                     
                     
                     
                  
                     
