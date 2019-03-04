library(shiny)
library(tidyverse)
library(readxl)
source('code/functions.R')

# Read in data
normalized_data_genelevel_tpm = read_csv("data/normalized_data_genelevel_tpm.csv")
geneSyns = read_tsv('data/GeneNames.tsv')

shinyServer(function(input, output) {
   
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
      ggplot(data, aes(Day, mean, colour = siRNA)) + geom_point(size = 4) + geom_path(size = 1) + 
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
      ggsave(file, width = plotDims()[1], height = plotDims()[2], units = c('in'))
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
  
})
