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
                  tabPanel("PSM", plotOutput("nPSMs"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize=10*1024^2)

  scp <- reactive({
    
    req(input$evidence_file)
    req(input$sample_annotation_file)
    req(input$label_suffix)
    
    evidence_data <- read.delim(input$evidence_file$datapath)
    meta_data <- read.delim(input$sample_annotation_file$datapath)
    
    scp_0 <- readSCP(featureData = evidence_data,
              colData = read.delim(input$sample_annotation_file$datapath),
              channelCol = "Channel",
              batchCol = "Raw.file",
              suffix = paste0(input$label_suffix, 1:nrow(meta_data)),
              removeEmptyCols = TRUE)
  
    # change zeros to NA, apply first filter
    scp_0 <- zeroIsNA(scp_0, 1:length(rowDataNames(scp_0)))
  
    req(input$PIF_cutoff)
    
    # apply first filter
    scp_0 <- filterFeatures(scp_0,
                            ~ Reverse != "+" &
                              Potential.contaminant != "+" &
                              !is.na(PIF) & PIF > PIF_cutoff())
  
  
    return(scp_0)
  })
  
  
  output$overview_plot <- renderPlot({
    if (!is.null(scp())) {
      plot(scp()) }
  })
  
  
    output$nPSMs <- renderPlot({
      scp_0 <- scp()
      print(scp_0)
      nPSMs <- dims(scp_0[1,])
      req(input$nPSM_bins)  
      ggplot(data.frame(nPSMs)) +
          aes(x = nPSMs) +
          geom_histogram(binwidth = nPSM_bins)
    })
}




# Run the application 
shinyApp(ui = ui, server = server)
