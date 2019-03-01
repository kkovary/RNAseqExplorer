library(shiny)
source('code/functions.R')

shinyUI(
  navbarPage('RNAseqExplorer',
             tabPanel('Time Course Plots',
                      sidebarPanel(
                        #img(src = 'logo.png', style = "float: left; width: 75px; margin-right: 10px; margin-top: 5px"),
                        titlePanel(strong("RNA-seq Time Course Plots")),
                        h6(em("Atefeh Rabiee, Nathan Abell, Mary N. Teruel")),
                        selectInput("plot_type", "Plot Type:", c("Mean","Loess")),
                        textOutput("result"),
                        
                        checkboxGroupInput("siRNA", "siRNA Data to Plot:",
                                           c("siCtrl" = "siNC","siPPARg" = "siPPARG"),
                                           selected = c("siNC","siPPARG"),
                                           inline = TRUE),
                            
                        checkboxGroupInput("time", "Time Points to Plot (Days):",
                                           c("0" = 0,"0.5" = 0.5, "1" = 1, "2" = 2, 
                                             "3" = 3, "4" = 4, "5" = 5, "6" = 6),
                                           selected = c(0,0.5,1,2,3,4,5,6),
                                           inline = TRUE),
                            
                        textInput("genes", "Genes (separate by space):", value = "Pparg"),
                            
                        actionButton('plot','Plot'),
                        
                        
                        h3(strong("Download")),
                            
                        # Download Plot Settings
                        h5("PDF Dimensions"),
                        splitLayout(
                          textInput("width", "Width (in)", value = 10),
                          textInput("height", "Height (in)", value = 5)
                          ),
                            
                        downloadButton("downloadPlot", "Export plot as PDF"),
                            
                        # Download Data Settings
                        downloadButton("downloadData", "Export table as CSV")
                            
                            
                        ),
                          
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Plots", 
                                   plotOutput("plot")
                                   ),
                          tabPanel("Data Table", 
                                   tableOutput("table")
                                   )
                        )
                        )
                      ),
             tabPanel('Volcano Plots')
             )
  )
