library(shiny)
library("scp")
library("ggplot2")
library("magrittr")
library("dplyr")

reactlog::reactlog_enable()

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Single Cell Proteomics Workbench"),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      # read in of the data
      fileInput("evidence_file", "Upload evidence.txt from MaxQuant", accept = c("text")),
      fileInput("sample_annotation_file", "Upload sample annotatio file", accept = c("text")),
      textInput("label_suffix", "Input the label suffix (e.g. _TMT"),
      
      # cutoff for PIF
      numericInput("PIF_cutoff", "Input cutoff value for parental ion fraction. PSMs larger than the value remain", 1, min = 0, max=1, step = 0.1),
      
      sliderInput("nPSM_bins", "select number of bins for nPSM Histogram", round = T, min = 1, max = 500, value = 250),
      numericInput("qvalue_cutoff", "Input cutoff value for q-value. PSMs with q-value smaller than the value remain", 0.01, min = 0, max=1, step = 0.1)
    ),
    
    # Main Panel for the output
    mainPanel(
      #create tabs
      tabsetPanel(type = "tabs",
                  tabPanel("Overview", plotOutput("overview_plot")),
                  tabPanel("PSM")
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize=10*1024^2)
  
  # handle the data
  evidence_data <- reactive({
    req(input$evidence_file)
    read.delim(input$evidence_file$datapath)
  })
  meta_data <- reactive({
    req(input$sample_annotation_file)
    read.delim(input$sample_annotation_file$datapath)
  })
  label <- reactive({
    req(input$label_suffix)
  })
  
  
  # read in the data into scp object
  scp <- reactive({
    readSCP(featureData = evidence_data(),
            colData = meta_data(),
            channelCol = "Channel",
            batchCol = "Raw.file",
            suffix = paste0(label(), 1:nrow(meta_data())),
            removeEmptyCols = TRUE)
  })
  
  
  output$overview_plot <- renderPlot({
    if (!is.null(scp())) {
      plot(scp()) }
  })
  
  # validation for variable
  PIF_cutoff <- reactive({
    req(input$PIF_cutoff)
  })
  
  # change zeros to NA, apply first filter
    scp <- reactive((zeroIsNA(scp(),
                 ~ Reverse != "+" &
                   Potential.contaminant != "+" &
                   !is.na(PIF) & PIF > PIF_cutoff())
        ))
  
  
  #  reactive({nPSMs <- dims(scp()[1,])})
  
  #  nPSM_bins <- reactive({
  #    req(input$nPSM_bins)
  #  })
  
  #  output$nPSMs <- renderPlot({
  #      ggplot(data.frame(nPSMs())) +
  #        aes(x = nPSMs()) +
  #        geom_histogram(binwidth = nPSM_bins())
  #  })
  
  
}




# Run the application 
shinyApp(ui = ui, server = server)