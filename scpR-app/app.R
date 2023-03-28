library(shiny)
library("scp")
library("ggplot2")
library("magrittr")
library("dplyr")
library("reshape2")
library(scater)
library(limma)


reactlog::reactlog_enable()

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Proteomics Workbench"),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      # the button to rule them all
      actionButton("update_button", "Press to run/update"),
      actionButton("help", "Help"),
      # read in of the data
      fileInput("evidence_file", "Upload evidence.txt from MaxQuant", accept = c("text")),
      fileInput("sample_annotation_file", "Upload sample annotation file", accept = c("text")),
      textInput("label_suffix", "Input the label suffix (e.g. _TMT)"),
      
      # cutoff for PIF
      numericInput("PIF_cutoff", "Input cutoff value for parental ion fraction. PSMs larger than the value remain", 0.1, min = 0, max=1, step = 0.1),
      
      numericInput("qvalue_cutoff", "Input cutoff value for q-value. PSMs with q-value smaller than the value remain", 0.01, min = 0, max=1, step = 0.0001),
      
      # minimal observation of a peptide within their corresponding razor protein
      numericInput("nObs_pep_razrpr", "Input minimal observation of a peptide within their corresponding razor peptide", 5, min = 1, max=100, step = 1),
      
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
                  tabPanel("Dimensionality reduction", plotOutput("PCA"), plotOutput("UMAP")),
                  tabPanel("Feature wise output", 
                           selectInput("selectedProtein", "Choose protein for observation", ""),
                           plotOutput("feature_subset")),
                  
                  tabPanel("Statistics",
                           actionButton("run_statistics", "Press to run statistics"),
                           actionButton("qqplot", "check dependencies for linear model"),
                           textOutput("model_vars"),
                           selectInput("model_design", "Choose your study design",
                                       choices = c("All pairwise comparison",
                                                   "Treatment versus Control",
                                                   "Additivity using a twofactor model",
                                                   "Interaction using a twofactor model")),
                           plotOutput("venn_diagram"),
                           plotOutput("volcano", hover = hoverOpts(id ="plot_hover")),
                           verbatimTextOutput("hover_info"),
                           tableOutput("protein_table")
                  )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize=10*1024^2)
  
  ## Data structures and objects
  # sample_annotation needs to be accessible for plots
  meta_data <- reactive({
    req(input$sample_annotation_file)
    meta_data_0 <- read.delim(input$sample_annotation_file$datapath)
    return(meta_data_0)
  })
  
  # result of analysis pipleline
  scp <- eventReactive(input$update_button, {
    withProgress(message= "running analysis", value=0, {
      req(input$evidence_file)
      req(meta_data)
      req(input$label_suffix)
      
      evidence_data <- read.delim(input$evidence_file$datapath)
      meta_data_0 <- meta_data()
      
      incProgress(1/17, detail = paste("read Data"))
      scp_0 <- readSCP(featureData = evidence_data,
                       colData = meta_data_0,
                       channelCol = "Channel",
                       batchCol = "Raw.file",
                       suffix = paste0(input$label_suffix, 1:nrow(meta_data_0)),
                       removeEmptyCols = TRUE)
      
      # # change zeros to NA, apply first filter
      incProgress(2/17, detail = paste("replacing zeros with NA"))
      scp_0 <- zeroIsNA(scp_0, 1:length(rowDataNames(scp_0)))
      
      req(input$PIF_cutoff)
      
      # apply first filter
      incProgress(3/17, detail = paste("filter contaminants and PIF"))
      scp_0 <- filterFeatures(scp_0,
                              ~ Reverse != "+" &
                                Potential.contaminant != "+" &
                                !is.na(PIF) & PIF > input$PIF_cutoff)
      
      
      # compute qvalues_PSMs to filter out by FDR
      incProgress(4/17, detail=paste("calculate q-value for PSMs"))
      scp_0 <- pep2qvalue(scp_0,
                          i = 1:length(rowDataNames(scp_0)),
                          PEP = "PEP", # by reference the dart_PEP value is used
                          rowDataName = "qvalue_PSMs")
      incProgress(5/17, detail=paste("calculate q value for proteins"))
      scp_0 <- pep2qvalue(scp_0,
                          i = 1:length(rowDataNames(scp_0)),
                          PEP = "PEP",
                          groupBy = "Leading.razor.protein",
                          rowDataName = "qvalue_proteins")
      
      req(input$qvalue_cutoff)
      incProgress(6/17, detail=paste("filter according to q-value"))
      scp_0 <- filterFeatures(scp_0, ~ qvalue_proteins < input$qvalue_cutoff)
      
      
      # aggregate PSMS to peptides
      incProgress(7/17, detail=paste("aggregating features"))
      scp_0 <- aggregateFeaturesOverAssays(scp_0,
                                           i = 1:length(rowDataNames(scp_0)),
                                           fcol = "Modified.sequence",
                                           name = paste0("peptides_", names(scp_0)),
                                           fun = matrixStats::colMedians, na.rm = TRUE)
      
      file_name <- meta_data_0$Raw.file[1]
      peptides <- paste("peptides_", as.character(file_name), sep = "")
      
      # calculate median reporter IO intensity
      incProgress(8/17, detail=paste("calculate reporter ion intensity"))
      medians <- colMedians(assay(scp_0[[peptides]]), na.rm = TRUE)
      scp_0$MedianRI <- medians
      
      # Filter based on the median CV -> remove covariant peptides over multiple proteins
      
      req(input$nObs_pep_razrpr)
      incProgress(9/17, detail=paste("calculate covariance per cell"))
      scp_0 <- medianCVperCell(scp_0,
                               i = length(rowDataNames(scp_0)),
                               groupBy = "Leading.razor.protein",
                               nobs = input$nObs_pep_razrpr, 
                               norm = "div.median",
                               na.rm = TRUE,
                               colDataName = "MedianCV")
      
      incProgress(10/17, detail=paste("filtering according to covariance"))
      req(input$MedCV_thresh)
      scp_0 <- scp_0[, !is.na(scp_0$MedianCV) & scp_0$MedianCV < input$MedCV_thresh, ]
      
      incProgress(11/17, detail=paste("normalizing peptide data"))
      # divide column by median
      scp_0 <-normalizeSCP(scp_0, 2, name="peptides_norm_col", method = "div.median")
      #divide rows my mean
      scp_0 <-normalizeSCP(scp_0, 3, name="peptides_norm", method = "div.mean")
      
      
      incProgress(12/17, detail=paste("remove peptides with high missing rate"))
      req(input$pNA)
      scp_0 <- filterNA(scp_0,
                        i = "peptides_norm",
                        pNA = input$pNA)
      
      
      incProgress(13/17, detail=paste("log-transforming peptide data")) 
      req(input$transform_base)
      scp_0 <- logTransform(scp_0,
                            base = as.integer(input$transform_base),
                            i = "peptides_norm",
                            name = "peptides_log")
      
      incProgress(14/17, detail=paste("aggregate peptides to proteins"))
      scp_0 <- aggregateFeatures(scp_0,
                                 i = "peptides_log",
                                 name = "proteins",
                                 fcol = "Leading.razor.protein",
                                 fun = matrixStats::colMedians, na.rm = TRUE)
      
      incProgress(15/17, detail=paste("normalizing proteins"))
      scp_0 <- sweep(scp_0, i = "proteins",
                     MARGIN = 2,
                     FUN = "-",
                     STATS = colMedians(assay(scp_0[["proteins"]]),
                                        na.rm = TRUE),
                     name = "proteins_norm_col")
      
      # Center rows with mean
      scp_0 <- sweep(scp_0, i = "proteins_norm_col",
                     MARGIN = 1,
                     FUN = "-",
                     STATS = rowMeans(assay(scp_0[["proteins_norm_col"]]),
                                      na.rm = TRUE),
                     name = "proteins_norm")
      
      
      sce <- getWithColData(scp_0, "proteins_norm")
      
      batch <- colData(sce)$Raw.file
      model <- model.matrix(~ SampleType, data = colData(sce))
      
      scp_0 <- addAssay(scp_0,
                        y = sce,
                        name = "proteins_dim_red")
      
      scp_0 <- addAssayLinkOneToOne(scp_0,
                                    from = "proteins_norm",
                                    to = "proteins_dim_red")
      
      
      
      
      incProgress(16/17, detail=paste("running dimensionality reduction"))
      
      scp_0[["proteins_dim_red"]] <- runPCA(scp_0[["proteins_dim_red"]],
                                            ncomponents = 5,
                                            ntop = Inf,
                                            scale = TRUE,
                                            exprs_values = 1,
                                            name = "PCA")
      
      scp_0[["proteins_dim_red"]] <- runUMAP(scp_0[["proteins_dim_red"]],
                                             ncomponents = 2,
                                             ntop = Inf,
                                             scale = TRUE,
                                             exprs_values = 1,
                                             n_neighbors = 3,
                                             dimred = "PCA",
                                             name = "UMAP")
      
      
      
      incProgress(17/17, detail=paste("analysis finish"))
    })
    return(scp_0)
  })
  
  # expression matrix is used by plot funcitons & stat pipeline
  exp_matrix <- reactive({
    scp_0 <- scp()
    exp_matrix_0 <- data.frame(assay(scp_0[["proteins_dim_red"]]))
    colnames(exp_matrix_0) <- scp_0$SampleType  
    return(exp_matrix_0)
  })
  
  # statistical pipeline
  stat_result <- eventReactive(input$run_statistics, {
    withProgress(message= "running statistical analysis", value=0, {
      incProgress(1/5, detail=paste("read data"))
      scp_0 <- scp()
      
      incProgress(2/5, detail=paste("reading expression matrix"))
      exp_matrix_0 <- exp_matrix()
      
      factorize_var <- function(i_vector) {
        fact_vect <- factor(i_vector)
        levels <- unique(fact_vect)
        num_values <- seq_along(levels)
        lookup_table <- data.frame(fact_vect = levels, num_var = num_values)
        
        num_var <- sapply(fact_vect, function(x) {
          lookup_table$num_var[lookup_table$fact_vect == x]
        })
        
        fact_vect <- factor(num_var)
        return(fact_vect)
      }
      
      factors_sample_type <- factorize_var(colnames(exp_matrix_0))
      
      
      incProgress(3/5, detail=paste("creating design matrix"))
      
      req(input$model_design)
      
      # Create a design matrix
      if (input$model_design == "All pairwise comparison") {
        # groupwise --> delete duplicates
        design <- model.matrix(~ 0+factors_sample_type)
        colnames(design) <- unique(scp_0$SampleType)
      }
      # treatment versus control
      else if (input$model_design == "Treatment versus Control") {
        design <- model.matrix(~factors_sample_type)
        colnames(design) <- c("Intercept", unique(scp_0$SampleType)[-1])
      }
      else if (input$model_design == "Additivity using a twofactor model") {
        factors_user_choice <- factorize_var(model_vals$factor_vector)
        design <- model.matrix(~factors_sample_type + factors_user_choice)
        user_colnames <- sprintf(paste(as.character(model_vals$vector_name),"[%s]"), seq(2:length(unique(factors_user_choice)))+1)
        colnames(design) <- c(unique(scp_0$SampleType), user_colnames)
      }
      else if (input$model_design == "Interaction using a twofactor model") {
        factors_user_choice <- factorize_var(model_vals$factor_vector)
        design <- model.matrix(~factors_sample_type * factors_user_choice)
        # user_colnames <- sprintf(paste(as.character(model_vals$vector_name),"[%s]"), seq(2:length(unique(factors_user_choice)))+1)
        # colnames(design) <- c(unique(scp_0$SampleType), user_colnames)
      }
      
      
      incProgress(4/5, detail=paste("Fit the expression matrix to a linear model"))
      fit <- lmFit(exp_matrix_0, design)
      incProgress(5/5, detail=paste("Bayes statistics of differential expression"))
      # *There are several options to tweak!*
      fit <- eBayes(fit)
    })
    return(fit)
  })
  
  # initial qq counter for the first plot 
  qq_count <- reactiveVal(1)
  
  ### graphical outputs and user interface
  
  ## reactive events and observers
  # Help Button and content of the help menu
  observeEvent(input$help,{
    showModal(modalDialog(easyClose = T, 
                          title = "Proteomics Workbench",
                          HTML("Created by Lukas Gamp.<br> <br>
      Workflow as described in the publication: <br>
      Multiplexed single-cell proteomics using SCoPE2 <br>
             10.1038/s41596-021-00616-z <br>
             <br>
             <br>
             1.) Read and filter Data <br>
             zeros are replaced with NA <br>
             Peptide sequence matches (PSMs) with parental ion fraction larger than the threshold will be kept. 
             Contaminants (reverse matches in the database) will be excluded. <br>
             PSMs will be shown in the summary barplot (first bar). <br>
             <br>
             <br>
             2.) Posterior error probability to q-value transformation and filtering <br>
             q-value reports significance of a match (against randomness). Adjust threshold accordingly <br>
             <br>
             <br>
             3.) PSMs are aggregated to peptides.<br>
             <br>
             <br>
             4.) Calculation of median reporter ion intensity (RI). <br>
             visualization in RI plot. <br>
             <br>
             <br>
             5.) Covariance across razor protein calculation and filtering according to n-observations (nobs) <br>
             Adapt nobs and accepted covariance according to your preference of analysis in order to filter out noisy quantification. Visualization can be found in the CV plot. <br>
             <br>
             <br>
             6.) Normalization of peptide data <br>
             Relative intensities will be divided by the median relative intensities. <br>
             <br>
             <br> 
             7.) Peptides with missing data deletion <br>
             Accepted threshold peptides with missing data can be adjusted as a fraction. <br>
             <br>
             <br>
             8.) Log transformation of the peptide data <br>
             Log base can be chosen (2 or 10) <br>
             <br>
             <br>
             9.) Peptide to protein aggregation <br>
             <br>
             <br>
             10.) Normalzation of protein data <br>
             Since log transformed --> subtractive method. <br>
             <br>
             <br>
             11.) Dimensionality reduction. <br>
             <br>
             <br>
             12.) Protein wise observation <br>
             select your protein of interest. <br>
             <br>
             <br>
             ")
                          
                          
                          
                          
    ))
  })
  
  # observer for protein list
  observe({
    updateSelectInput(session, "selectedProtein", choices = protein_list())
  })
  
  # reactive element for the protein list --> update if the scp object changes
  protein_list <- reactive({
    scp_0 <- scp()
    list <- rowData(scp_0)[["proteins"]][,1]
    return(list)
  })
  
  # reactive element for hover function of volcano plot
  displayed_text <- reactive({
    req(input$plot_hover)
    hover <- input$plot_hover
    dist <- sqrt((hover$x - stat_result()$coef)^2 + (hover$y - -log10(stat_result()$p.value))^2)
    
    if(min(dist) < 0.3) {
      rownames(stat_result())[which.min(dist)]
    } else {
      NULL
    }
  })

  # observer for qqplot interface
  observeEvent(input$qqplot, {
    showModal(qqModal())
  })
  
  
  # observer for the next button to increment the counter and show next qq
  observeEvent(input$next_plot, {
    exp_matrix_0 <- exp_matrix()
    sample_types <- unique(colnames(exp_matrix_0))
    
    i <- qq_count()
    
    if(i < length(sample_types)) {
      j <- i + 1
    }
    else {
      j <- 1
    }
    
    qq_count(j)
  })  
  
  # observer for the previous button to increment the counter and show previous qq
  observeEvent(input$previous_plot, {
    exp_matrix_0 <- exp_matrix()
    sample_types <- unique(colnames(exp_matrix_0))
    
    i <- qq_count()
    
    if(i > 1) {
      j <- i - 1
    }
    else {
      j <- 1
    }
    
    qq_count(j)
  })  
  
  ## Plots 
  # pathway of the data first tab
  output$overview_plot <- renderPlot({
    if (!is.null(scp())) {
      plot(scp()) }
  })
  
  # summary barchart second tab
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
  
  # Boxplot of reporter ion intensity third tab
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
  
  # covariance across razor peptides fourth tab
  output$CV_median <- renderPlot({
    req(scp())
    req(meta_data())
    
    scp_0 <- scp()
    meta_data_0 <- meta_data()
    
    file_name <- meta_data_0$Raw.file[1]
    peptides <- paste("peptides_", as.character(file_name), sep = "")
    
    getWithColData(scp_0, peptides) %>%
      colData %>%
      data.frame %>%
      ggplot(aes(x = MedianCV,
                 fill = SampleType)) +
      geom_boxplot()
  })
  
  # principle component analysis in fifth tab
  output$PCA <- renderPlot({
    scp_0 <- scp()
    plotReducedDim(scp_0[["proteins_dim_red"]],
                   dimred = "PCA",
                   colour_by = "SampleType",
                   point_alpha = 1,
                   point_size=3)
  })
  
  # Umap dimensionality reduction in fith tab
  output$UMAP <- renderPlot({
    scp_0 <- scp()
    plotReducedDim(scp_0[["proteins_dim_red"]],
                   dimred = "UMAP",
                   colour_by = "SampleType",
                   point_alpha = 1,
                   point_size = 3)
  })
  
  # protein wise visualisation 
  output$feature_subset <- renderPlot({
    scp_0 <- scp()
    label_suffix <- input$label_suffix
    subsetByFeature(scp_0, input$selectedProtein) %>%
      ## Format the `QFeatures` to a long format table
      longFormat(colvars = c("Raw.file", "SampleType", "Channel")) %>%
      data.frame %>%
      ## This is used to preserve ordering of the samples and assays in ggplot2
      mutate(assay = factor(assay, levels = names(scp_0)),
             Channel = sub("Reporter.intensity.", label_suffix, Channel),
             Channel = factor(Channel, levels = unique(Channel))) %>%
      ## Start plotting
      ggplot(aes(x = Channel, y = value, group = rowname, col = SampleType)) +
      geom_point() +
      ## Plot every assay in a separate facet
      facet_wrap(facets = vars(assay), scales = "free_y", ncol = 3) +
      ## Annotate plot
      xlab("Channels") +
      ylab("Intensity (arbitrary units)") +
      ## Improve plot aspect
      theme(axis.text.x = element_text(angle = 90),
            strip.text = element_text(hjust = 0),
            legend.position = "bottom")
  })
  
  # volcanoplot in statistic tab
  output$volcano <- renderPlot({
    volcanoplot(stat_result(), cex = 0.5) 
  })
  
  # hover output below volcano chart
  output$hover_info <- renderPrint({
    req(displayed_text())
    
    cat("UniProt-ID\n")
    displayed_text()
  })
  
  # cutoff for the p-adj, foldchange
  # qlucor 
  # contrasts 
  
  # show significant proteins in the table
  output$protein_table <- renderTable({
    topTable(stat_result(), number = 10, adjust = "BH")
  }, rownames = T)
  
  # venn diagram for significant proteins
  output$venn_diagram <- renderPlot({
    req(stat_result())
    vennDiagram(decideTests(stat_result()))
  })
  
  output$model_vars <- renderText({
    req(scp())
    unique(scp()$SampleType)
  })
  
  model_vals <- reactiveValues(data = NULL)

  dataModal <- function(failed_1 = FALSE, failed_2=FALSE) {
    modalDialog(
      textInput("factor_vector", "Insert character or number vector in the length of your sample size (comma delim.)",
                placeholder = 'Try "1,2,3,4,5,6,1,2,3..." or "C,C,T,T..."'
      ),
      span('(Try a valid vector like "1,2,3,4", ',
           'then a non-valid vector like "abc-cde")'),
      if (failed_1)
        div(tags$b("Invalid vector!", style = "color: red;")),
      
      textInput("vector_name", "Insert name of your vector"),
      if (failed_2)
        div(tags$b("Name your vector!", style = "color: red;")),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ok", "OK") 
      )
    )
  }
  
  # plot the qq_count-th column
  output$qqPlot <- renderPlot({
    exp_matrix_0 <- exp_matrix()
    sample_types <- unique(colnames(exp_matrix_0))
    i <- qq_count()
    qqnorm(exp_matrix_0[[sample_types[i]]], main=paste("QQ-plot of ", sample_types[i]))
    qqline(exp_matrix_0[[sample_types[i]]])
  })
  
  output$hist <- renderPlot({
    exp_matrix_0 <- exp_matrix()
    sample_types <- unique(colnames(exp_matrix_0))
    i <- qq_count()
    hist(exp_matrix_0[[sample_types[i]]], main=paste("Histogram of ", sample_types[i]), xlab= "")
  })
  
  # create interface for the qqmodal dialog
  qqModal <- function() {
    modalDialog(
      actionButton("previous_plot", "Show previous"),
      actionButton("next_plot", "Show next"),
      plotOutput("qqPlot"),
      plotOutput("hist"),
      footer = tagList(
        modalButton("Dismiss")
      )
    )
  }
  
  observeEvent(input$model_design, {
    if (input$model_design == "Additivity using a twofactor model") {
      showModal(dataModal())
    }
    else if (input$model_design == "Interaction using a twofactor model") {
      showModal(dataModal())
    }
  })
  
  observeEvent(input$ok, {
    factor_vect <- strsplit(input$factor_vector, ",\\s*")[[1]]
    
    if (!is.null(factor_vect) 
        && is.vector(factor_vect) 
        && length(factor_vect) == length(scp()$SampleType)) {
      model_vals$factor_vector <- factor_vect
      removeModal()
    } else {
      showModal(dataModal(failed_1 = TRUE))
    }
    
    if (!is.null(input$vector_name)) {
      model_vals$vector_name <- input$vector_name
      removeModal()
    } else {
      showModal(dataModal(failed_2 = TRUE))
    }
  })
  
  
}




# Run the application 
shinyApp(ui = ui, server = server)