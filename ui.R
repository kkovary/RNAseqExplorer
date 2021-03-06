library(shiny)
library(readr)
source('code/functions.R')

volData = read_csv('data/volData.csv')

shinyUI(
  navbarPage('RNAseqExplorer',
             tabPanel('Time Course Plots',
                      sidebarPanel(
                        #img(src = 'logo.png', style = "float: left; width: 75px; margin-right: 10px; margin-top: 5px"),
                        titlePanel(strong("RNA-seq Time Course Plots")),
                        h6(em("Atefeh Rabiee, Nathan Abell, Mary N. Teruel")),
                        selectInput("plot_type", "Plot Type:", c("Mean","Loess")),
                        sliderInput('pointSize', 'Point Size', min = 0, max = 5, value = 2.5, step = 0.1),
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
             tabPanel('Volcano Plots',
                      sidebarPanel(
                        h6(strong('Please be patient,'), em('it may take a minute to display the plot.')),
                        tags$hr(),
                        tags$p("To create volcano plots, select two conditions to compare.
                               For example, to see how genes are differentially expressed
                               from day 0 to day 6, select NC_6 as the numerator, and NC_0
                               as the denominator (NC is non-targeting siRNA, PPARG is 
                               Pparg targeting siRNA."),
                        selectInput("volNumerator", "Numerator:", choices = c(NA,unique(volData$Condition)), selected = 'NC_6'),
                        selectInput("volDenominator", "Denominator:", choices = c(NA,unique(volData$Condition)), selected = 'NC_0'),
                        tags$hr(),
                        tags$p('Select the pvalue and fold change cutoff'),
                        numericInput('pvalCut', 'pvalue Cutoff', value = 0.05, min = 0, max = 1, step = 0.01),
                        numericInput('fcCut', 'Fold Change Cutoff', value = 2, step = 0.25),
                        actionButton('volPlotButton','Plot'),
                        tags$hr(),
                        tags$p(strong('Download')),
                        h5("PDF Dimensions"),
                        splitLayout(
                          textInput("volWidth", "Width (in)", value = 10),
                          textInput("volHeight", "Height (in)", value = 5)
                        ),
                        
                        downloadButton("volDownloadPlot", "Export plot as PDF"),
                        
                        # Download Data Settings
                        downloadButton("volDownloadData", "Export table as CSV")
                        
                      ),
                      mainPanel(
                        plotOutput("volPlot"),
                        verbatimTextOutput("head")
                      )
             ),
             tabPanel('Dynamics Categories',
                      sidebarPanel(
                        checkboxGroupInput('d0_siRNA', 
                                           label = 'Day 0 siCtrl/siPparg',
                                           choices = c('+','=','-'),
                                           inline = TRUE),
                        checkboxGroupInput('dm_siRNA', 
                                           label = 'Dex+IBMX siCtrl/siPparg',
                                           choices = c('+','=','-'),
                                           inline = TRUE),
                        checkboxGroupInput('dm_NC', 
                                           label = 'Dex+IBMX siCtrl DM/0',
                                           choices = c('+','=','-'),
                                           inline = TRUE),
                        checkboxGroupInput('dm_PPARG', 
                                           label = 'Dex+IBMX siPparg DM/0',
                                           choices = c('+','=','-'),
                                           inline = TRUE),
                        checkboxGroupInput('ins_siRNA', 
                                           label = 'Ins siCtrl/siPparg',
                                           choices = c('+','=','-'),
                                           inline = TRUE),
                        checkboxGroupInput('ins_NC', 
                                           label = 'Ins siCtrl DM/0',
                                           choices = c('+','=','-'),
                                           inline = TRUE),
                        checkboxGroupInput('ins_PPARG', 
                                           label = 'Ins siPparg DM/0',
                                           choices = c('+','=','-'),
                                           inline = TRUE),
                        
                        actionButton('dynPlotButton','Plot')
                      ),
                      mainPanel(
                        plotOutput('dynPlot', height = '2000px'),
                        #verbatimTextOutput("head")
                      )
             )
  )
)