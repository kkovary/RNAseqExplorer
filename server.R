library(shiny)
library(tidyverse)
library(readxl)
library(httr)
library(viridis)
library(scales)
library(gtools)
source('code/functions.R')

# Read in data
normalized_data_genelevel_tpm = read_csv("data/normalized_data_genelevel_tpm.csv")
geneSyns = read_csv('data/GeneNames.csv')
volData = read_csv('data/volData.csv')
norm <- read_csv('data/norm.csv')
cat_quant_long <- read_csv('data/cat_quant_long.csv')

shinyServer(function(input, output, session) {
  
  # Reads genes and formats them when plot buttion is pushed
  geneNames <- reactive({
    input$plot
    isolate(input$genes %>% strsplit(' ') %>% unlist() %>% tolower())
  })
  
  # Formats the data for plotting
  datasetInput <- reactive({
    geneSearch(geneNames(), geneSyns, normalized_data_genelevel_tpm) %>% 
      filter(siRNA %in% input$siRNA, Day %in% input$time)
  })
  
  # Formats the data for table output
  tableFormat <- reactive({
    
    tab = datasetInput() %>% unite("Day_Replicate", c("Day", "Replicate")) %>% 
      spread("Day_Replicate","TPM") %>% unite("Condition", c("GeneName","siRNA"))
    
    mat = t(tab[2:ncol(tab)])
    mat = cbind(colnames(tab[,2:ncol(tab)]),mat)
    colnames(mat) = c("Sample",as.character(tab$Condition))
    as.data.frame(mat) %>% separate(Sample, into = c("Day","Replicate"), sep = "\\_")
  })
  
  # Updates dropdown menu with searched genes
  observe({
    updateSelectInput(session = session, inputId = "uniprot", choices = geneNames())
  })
  
  # Write Uniprot HTML File
  oijfewoji <- reactive({
    selected = filter(geneSyns, tolower(GN_Syn) == tolower(input$uniprot))
    selected = selected[1,1] %>% as.character()
    
    paste0('https://www.uniprot.org/uniprot/',selected)
  })
  
  output$Uniprot <- renderText({
    #request <- GET(url = url())
    #writeLines(content(request, as="text"), file('uniprot.html'))
    #content(request, as = "text")
    as.character(url())
  })
  
  
  # Defines the different plots
  datasetPlot <- reactive({
    if(input$plot_type == 'Loess'){
      ggplot(datasetInput(), aes(Day, TPM, colour = siRNA)) + geom_point() + geom_smooth(method = loess) +
        facet_wrap(~GeneName, scales = "free") + scale_x_continuous(breaks = as.numeric(input$time)) + 
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_manual(values=c("#ED1C24", "#213F99"))
    }
    
    else{
      data = datasetInput() %>% group_by(GeneName,siRNA,Day) %>% summarise(mean = mean(TPM, na.rm = T), stdev = sd(TPM, na.rm = T))
      ggplot(data, aes(Day, mean, colour = siRNA)) + geom_point(size = input$pointSize) + geom_path(size = 1) + 
        geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev), width = 0.25, size = 1) +
        facet_wrap(~GeneName, scales = "free") + scale_x_continuous(breaks = as.numeric(input$time)) + 
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_manual(values=c("#ED1C24", "#213F99")) + ylab('Mean TPM')
    }
  })
  
  plotDims <- reactive({
    c(as.numeric(input$width), as.numeric(input$height))
  })
  
  # Plot Data  
  output$plot <- renderPlot({
    
    datasetPlot()
    
  })
  
  # Download PDF of plotted dataset ----
  output$downloadPlot <-downloadHandler(
    filename = function() {
      "Plot.pdf"
    },
    content = function(file){
      #x = plotDims()
      renderPlot({
        datasetPlot()
      })
      ggsave(file, width = plotDims()[1], height = plotDims()[2], units = c('in'), useDingbats=FALSE)
    }
  )
  
  # Table of selected dataset ----
  output$table <- renderTable({
    tableFormat()
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$genes, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(tableFormat(), file, row.names = FALSE)
    }
  )
  
  # Volcano Plot
  volPlotData <- reactive({
    volPlotDataFun(input$volNumerator, input$volDenominator, volData, input$pvalCut, input$fcCut) %>%
      filter(!is.na(Hit))
  })
  
  output$volPlot <- renderPlot({
    ggplot(volPlotData(), aes(x = log2(foldChange), y = -log10(pvalue), colour = Hit)) + 
      geom_vline(xintercept = log2(c(input$fcCut, 1/input$fcCut)), colour = 'red', linetype = 'dashed') + 
      geom_hline(yintercept = -log10(input$pvalCut), colour = 'red', linetype = 'dashed') + 
      geom_point(alpha = 0.25) + theme_bw() + scale_color_manual(values = c('#bababa','#e08214')) +
      theme(legend.position="none")
  })
  output$head <- renderPrint(
    volPlotData() %>% head()
  )
  # Download PDF of Volcano Plot
  output$volDownloadPlot <-downloadHandler(
    filename = function() {
      paste(input$volNumerator,' vs ',input$volDenominator, ".pdf", sep = "")
    },
    content = function(file){
      #x = plotDims()
      renderPlot({
        volPlot()
      })
      ggsave(file, width = plotDims()[1], height = plotDims()[2], units = c('in'))
    }
  )
  
  # Downloadable csv Volcano Plot
  output$volDownloadData <- downloadHandler(
    filename = function() {
      paste(input$volNumerator,' vs ',input$volDenominator, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(volPlotData(), file, row.names = FALSE)
    }
  )
  
  # Dynamics heatmap
  dynPlot <- eventReactive(
    input$dynPlotButton, 
    {
      cat_heatmap( 
        d0_siRNA_ = input$d0_siRNA, 
        dm_siRNA_ = input$dm_siRNA, 
        dm_NC_ = input$dm_NC,
        dm_PPARG_ = input$dm_PPARG,
        ins_siRNA_ = input$ins_siRNA,
        ins_NC_ = input$ins_NC,
        ins_PPARG_ = input$ins_PPARG,
        output = 'heatmap'
      )
      
    }
  )
  
  output$dynPlot <- renderPlot({
    dynPlot()
  })
  
})
