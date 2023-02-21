library(shiny)
library("scp")
library("ggplot2")
library("magrittr")
library("dplyr")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Single Cell Proteomics Workbench"),

    # Sidebar
    sidebarLayout(
        sidebarPanel(
            fileInput("evidence_file", "Upload evidence.txt from MaxQuant",
                      accept = c("text")),
            fileInput("sample_annotation_file", "Upload sample annotatio file",
                      accept = c("text")),
            textInput("label_suffix", "Input the label suffix (e.g. _TMT"),
            numericInput("PIF_cutoff", "Input cutoff value for parental ion fraction. PSMs larger than the value remain", 1, min = 0, max=1, step = 0.1),
            sliderInput("nPSM_bins", "select number of bins for nPSM Histogram", round = T, min = 1, max = 500, value = 250),
            numericInput("qvalue_cutoff", "Input cutoff value for q-value. PSMs with q-value smaller than the value remain", 0.01, min = 0, max=1, step = 0.1)
            
        ),
        
        # Main Panel for the output
        mainPanel(
          
          #create tabs
          
          tabsetPanel(type = "tabs",
                      tabPanel("Overview", plotOutput("overview")),
                      tabPanel("PSM", plotOutput("nPSMs"))
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  observe({
    # input declaration
    evidence_file = input$evidence_file
    sample_annotation_file = input$sample_annotation_file
    label_suffix = input$label_suffix
    PIF_cutoff = input$PIF_cutoff
    qvalue_cutoff = input$qvalue_cutoff
    # check if null
    if (is.null(evidence_file) || is.null(sample_annotation_file) || is.null(label_suffix) || is.null(PIF_cutoff) || is.null(qvalue_cutoff)){
      return(NULL)
    }
    
    # handle the input
    mqScpData = read.delim(evidence_file$datapath)
    sampleAnnotation = read.delim(sample_annotation_file$datapath)
    
    # read in the data into scp object
    scp <- readSCP(featureData = mqScpData,
                   colData = sampleAnnotation,
                   channelCol = "Channel",
                   batchCol = "Raw.file",
                   suffix = paste0(label_suffix, 1:nrow(sampleAnnotation)),
                   removeEmptyCols = TRUE)
    
    # convert Zeros to NA 
    scp <- zeroIsNA(scp, i=1:length(rowDataNames(scp)))
    
    scp <- filterFeatures(scp,
                          ~ Reverse != "+" &
                            Potential.contaminant != "+" &
                            !is.na(PIF) & PIF > PIF_cutoff)
    
    
    
    
    output$overview <- renderPlot({
      plot(scp)
    })
    
  # increase performance and check only necessary values within observe
    observe({
      
      # input declaration
      
      nPSM_bins = input$nPSM_bins
      
      # check if null
      if (is.null(nPSM_bins)){
        return(NULL)
      }
      
      nPSMs <- dims(scp)[1, ]
      output$nPSMs <- renderPlot({
      ggplot(data.frame(nPSMs)) +
        aes(x = nPSMs) +
        geom_histogram(binwidth = nPSM_bins)
      })
    })
    
  })
}




# Run the application 
shinyApp(ui = ui, server = server)
