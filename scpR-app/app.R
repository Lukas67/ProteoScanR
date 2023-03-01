library(shiny)
library("scp")
library("ggplot2")
library("magrittr")
library("dplyr")
library("reshape2")


reactlog::reactlog_enable()

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Proteomics Workbench"),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      # the button to rule them all
      actionButton("update_button", "Press to update"),
      
      # read in of the data
      fileInput("evidence_file", "Upload evidence.txt from MaxQuant", accept = c("text")),
      fileInput("sample_annotation_file", "Upload sample annotation file", accept = c("text")),
      textInput("label_suffix", "Input the label suffix (e.g. _TMT"),
      
      # cutoff for PIF
      numericInput("PIF_cutoff", "Input cutoff value for parental ion fraction. PSMs larger than the value remain", 0.1, min = 0, max=1, step = 0.1),
      
      numericInput("qvalue_cutoff", "Input cutoff value for q-value. PSMs with q-value smaller than the value remain", 0.01, min = 0, max=1, step = 0.0001),
      
      # minimal observation of a peptide within their corresponding razor protein
      numericInput("nObs_pep_razrpr", "Input minimal observation of a peptide within their corresponding razor protein", 5, min = 1, max=100, step = 1),
      
      # max covariance accepted
      numericInput("MedCV_thresh", "Input maximal Covariance accepted", 0.65 , min = 0, max=1, step = 0.01),

      # remove peptides with missing data
      numericInput("pNA", "Input percentage threshold for peptide data", 0.99 , min = 0, max=1, step = 0.01),
      
      # choose log transform base
      selectInput("transform_base", "Choose log base for peptide data transformation", choices = c("2", "10"))
    ),
    
    # Main Panel for the output
    mainPanel(
      #create tabs
      tabsetPanel(type = "tabs",
                  tabPanel("Overview", plotOutput("overview_plot")),
                  tabPanel("Summary Barplot", plotOutput("summary_bar")),
                  tabPanel("Reporter Ion Intensity", plotOutput("RI_intensity")),
                  tabPanel("Covariance across razor peptides", plotOutput("CV_median")),
                  tabPanel("Missing Proteins", plotOutput("missing"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize=10*1024^2)

  scp <- eventReactive(input$update_button, {
    
    req(input$evidence_file)
    req(input$sample_annotation_file)
    req(input$label_suffix)
    
    evidence_data <- read.delim(input$evidence_file$datapath)
    meta_data <- read.delim(input$sample_annotation_file$datapath)
    
    scp_0 <- readSCP(featureData = evidence_data,
              colData = meta_data,
              channelCol = "Channel",
              batchCol = "Raw.file",
              suffix = paste0(input$label_suffix, 1:nrow(meta_data)),
              removeEmptyCols = TRUE)
  
    # # change zeros to NA, apply first filter
    scp_0 <- zeroIsNA(scp_0, 1:length(rowDataNames(scp_0)))

    req(input$PIF_cutoff)
   
    # apply first filter
    scp_0 <- filterFeatures(scp_0,
                            ~ Reverse != "+" &
                              Potential.contaminant != "+" &
                              !is.na(PIF) & PIF > input$PIF_cutoff)

    
    # compute qvalues_PSMs to filter out by FDR
    
    scp_0 <- pep2qvalue(scp_0,
                      i = 1:length(rowDataNames(scp_0)),
                      PEP = "PEP", # by reference the dart_PEP value is used
                      rowDataName = "qvalue_PSMs")
    
    scp_0 <- pep2qvalue(scp_0,
                      i = 1:length(rowDataNames(scp_0)),
                      PEP = "PEP",
                      groupBy = "Leading.razor.protein",
                      rowDataName = "qvalue_proteins")
    
    req(input$qvalue_cutoff)
    scp_0 <- filterFeatures(scp_0, ~ qvalue_proteins < input$qvalue_cutoff)
    
    
    # aggregate PSMS to peptides
    
    scp_0 <- aggregateFeaturesOverAssays(scp_0,
                                       i = 1:length(rowDataNames(scp_0)),
                                       fcol = "Modified.sequence",
                                       name = paste0("peptides_", names(scp_0)),
                                       fun = matrixStats::colMedians, na.rm = TRUE)
    
    file_name <- meta_data$Raw.file[1]
    peptides <- paste("peptides_", as.character(file_name), sep = "")
    
    # calculate median reporter IO intensity
    
    medians <- colMedians(assay(scp_0[[peptides]]), na.rm = TRUE)
    scp_0$MedianRI <- medians
    
    # Filter based on the median CV -> remove covariant peptides over multiple proteins
    
    req(input$nObs_pep_razrpr)
    scp_0 <- medianCVperCell(scp_0,
                           i = length(rowDataNames(scp_0)),
                           groupBy = "Leading.razor.protein",
                           nobs = input$nObs_pep_razrpr, 
                           norm = "div.median",
                           na.rm = TRUE,
                           colDataName = "MedianCV")
    
    req(input$MedCV_thresh)
    scp_0 <- scp_0[, !is.na(scp_0$MedianCV) & scp_0$MedianCV < input$MedCV_thresh, ]
    
    # divide column by median
    scp_0 <-normalizeSCP(scp_0, 2, name="peptides_norm_col", method = "div.median")
    #divide rows my mean
    scp_0 <-normalizeSCP(scp_0, 3, name="peptides_norm", method = "div.mean")
    
    # remove peptides with high missing rate
    req(input$pNA)
    scp_0 <- filterNA(scp_0,
                    i = "peptides_norm",
                    pNA = input$pNA)
    
    
    # log-transform 
    req(input$transform_base)
    scp_0 <- logTransform(scp_0,
                        base = as.integer(input$transform_base),
                        i = "peptides_norm",
                        name = "peptides_log")
    
    # aggregate peptides to proteins
    scp_0 <- aggregateFeatures(scp_0,
                             i = "peptides_log",
                             name = "proteins",
                             fcol = "Leading.razor.protein",
                             fun = matrixStats::colMedians, na.rm = TRUE)
    
    # center cols
    scp_0 <- sweep(scp_0, i = "proteins",
                 MARGIN = 2,
                 FUN = "-",
                 STATS = colMedians(assay(scp_0[["proteins"]]),
                                    na.rm = TRUE),
                 name = "proteins_norm_col")
    
    ## Center rows with mean
    scp_0 <- sweep(scp_0, i = "proteins_norm_col",
                 MARGIN = 1,
                 FUN = "-",
                 STATS = rowMeans(assay(scp_0[["proteins_norm_col"]]),
                                  na.rm = TRUE),
                 name = "proteins_norm")
    
    # # missing value imputation
    # # bioc imputes feature space wise
    # scp_0 <- impute(scp_0,
    #               i = "proteins_norm",
    #               name = "proteins_imptd",
    #               method = "knn",
    #               k = 3, rowmax = 1, colmax= 1,
    #               maxp = Inf, rng.seed = 1234)
    # 
    # # can also be done with the Scope2 variant of KNN
    # # imputes by sample space 
    # scp <- imputeKnnSCoPE2(scp,
    #                        i = "proteins_norm2",
    #                        name = "proteins_impd",
    #                        k = 3)
    
    return(scp_0)
  })
  
  
  output$overview_plot <- renderPlot({
    if (!is.null(scp())) {
      plot(scp()) }
  })
  

    output$summary_bar <- renderPlot({
      scp_0 <- scp()
      
      nPSMs <- dims(scp_0)[1, 1]
      nPeps <- dims(scp_0)[1, 2]
      nProts <- dims(scp_0)[1, 3]
      summary_count <- data.frame(nPSMs, nPeps, nProts)
      summary_count <- reshape2::melt(summary_count)

      ggplot(summary_count, aes(x = variable, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        labs(x = "Number of PSMs, Peptides and Proteins", y = "Counts") +
        scale_fill_manual(values = c("darkblue", "darkgreen", "darkred"))
    })
    
    output$RI_intensity <- renderPlot({
      scp_0 <- scp()
      
      colData(scp_0) %>%
      data.frame %>%
      ggplot() +
      aes(x = MedianRI, 
          y = SampleType,
          fill = SampleType) +
      geom_boxplot() +
      scale_x_log10()
    })
    
    output$CV_median <- renderPlot({
      scp_0 <- scp()
      getWithColData(scp_0, peptides) %>%
        colData %>%
        data.frame %>%
        ggplot(aes(x = MedianCV,
                   fill = SampleType)) +
        geom_boxplot()
    })
    
    output$missing <- renderPlot({
      scp_0 <- scp()
      longFormat(scp_0[, , "proteins_norm"]) %>%
        data.frame %>%
        group_by(colname) %>%
        summarize(missingness = mean(is.na(value))) %>%
        ggplot(aes(x = missingness)) +
        geom_histogram()
      
    })
    
    
}




# Run the application 
shinyApp(ui = ui, server = server)
