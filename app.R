library(shiny)
library(tidyverse)

# Read in data
normalized_data_genelevel_tpm = read_csv("data/normalized_data_genelevel_tpm.csv")

# Define UI ----
ui <- fluidPage(
  titlePanel(strong("RNA-seq Time Course Plots")),
  h6(em("Atefeh Rabiee, Nathan Abell, Mary N. Teruel")),
  sidebarLayout(
    sidebarPanel(
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
      plotOutput("plot"),
      tableOutput("table")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  # Reactive value for selected dataset ----
  datasetInput <- reactive({
    plot_data <- filter(normalized_data_genelevel_tpm, tolower(GeneName) %in% tolower(unlist(strsplit(input$genes, " "))))
    plot_data <- gather(plot_data, "Sample", "TPM", 2:ncol(plot_data))
    plot_data <- plot_data %>% separate(Sample, into = c("siRNA", "Day", "Replicate"), sep = "\\_")
    plot_data$siRNA <- paste0("si",plot_data$siRNA)
    plot_data$Day <- as.numeric(plot_data$Day)
    plot_data <- filter(plot_data, siRNA %in% input$siRNA, Day %in% input$time)
  })
  
  tableFormat <- reactive({
    tab = datasetInput() %>% unite("Day_Replicate", c("Day", "Replicate")) %>% 
      spread("Day_Replicate","TPM") %>% unite("Condition", c("GeneName","siRNA"))
    
    mat = t(tab[2:ncol(tab)])
    mat = cbind(colnames(tab[,2:ncol(tab)]),mat)
    colnames(mat) = c("Sample",as.character(tab$Condition))
    as.data.frame(mat) %>% separate(Sample, into = c("Day","Replicate"), sep = "\\_")
  })
  

  
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
}

# Run the app ----
shinyApp(ui = ui, server = server)
