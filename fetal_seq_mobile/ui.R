

######################################
###                                ###
###         ui.R file for          ###
###       human-fetal-seq          ###
###          mobile                ###
###        application             ###
###                                ###
######################################

ui <- f7Page(title="Human Fetal-Seq Mobile",
             init=f7Init(theme = "light", skin="ios"),
             f7TabLayout(
               panels=tagList(
                 f7Panel(side="left", theme="dark", effect="reveal",
                         ##f7SmartSelect('entered_genes', label = 'Search for your gene of interest', choices=NULL, type=c("popup"), smart=TRUE, multiple=FALSE))),
                         selectizeInput('entered_genes', multiple = TRUE, choices = NULL, label = 'Search for gene of interest',
                                        options = list(maxItems = 1, placeholder = 'Example: Cx3cr1')))),
               navbar=f7Navbar(
                 title="Human Fetal-Seq Mobile", left_panel=TRUE, right_panel=FALSE),
               f7Tabs(
                 animated=TRUE, id='tabs',
                 f7Tab(tabName="Landing", active=TRUE, icon=f7Icon('home_fill'),
                       column(width=12,
                              tags$style("#landing_text{color: black; font-size: 20px;}"),
                              htmlOutput('landing_text'))),
                 f7Tab(tabName="Background", active=TRUE, icon = f7Icon("email"),
                       column(width=12,
                              htmlOutput('background_text'))),
                 f7Tab(tabName="Correlation Plots", active=TRUE, icon=f7Icon('layers'),
                       column(width=4,
                              withSpinner(plotlyOutput('plot1')))),
                 f7Tab(tabName="Sex Differences", active=TRUE, icon=f7Icon('graph_square', lib='ios'),
                       column(width=4,
                              withSpinner(plotlyOutput('sex_diff_plot')))),
                 f7Tab(tabName="Data Table", active=TRUE, icon=f7Icon('layers'),
                       column(width=4,
                              withSpinner(plotlyOutput('table1'))))
                       
               )
               )
               )





                     
                     
                     
                  
                     
