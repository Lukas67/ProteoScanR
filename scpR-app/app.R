library("shiny")
library("shinyWidgets")
library("shinydashboard")

library("scp")
library("SingleCellExperiment")
library("ggplot2")
library("magrittr")
library("dplyr")
library("reshape2")
library("scater")
library("limma")
library("CONSTANd")
library("stats")
library("impute")
library("sva")
library("plotly")

reactlog::reactlog_enable()

# Define UI for application that draws a histogram
ui <- fluidPage(
  # useShinyjs(),
  # Application title
  titlePanel("Proteomics Workbench"),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      # the button to rule them all
      actionButton("update_button", "Press to run/update"),
      actionButton("help", "Help"),
      
      # read in of the data
      materialSwitch(inputId = "file_level", value = F, label="Handle preprocessed expression set", status = "danger"),
      fileInput("evidence_file", "Upload evidence.txt from MaxQuant", accept = c("text")),
      fileInput("sample_annotation_file", "Upload sample annotation file", accept = c("text")),
      
      conditionalPanel(condition = "!input.file_level", 
                       selectInput("selectedSampleType_to_exclude", "Select SampleType to exclude", "", multiple=T),
                       # cutoff for PIF
                       numericInput("PIF_cutoff", "Input cutoff value for parental ion fraction. PSMs larger than the value remain", 0.1, min = 0, max=1, step = 0.1),
                       
                       numericInput("qvalue_cutoff", "Input cutoff value for q-value. PSMs with q-value smaller than the value remain", 0.01, min = 0, max=1, step = 0.0001),
                       
                       # minimal observation of a peptide within their corresponding razor protein
                       numericInput("nObs_pep_razrpr", "Input minimal observation of a peptide within their corresponding razor peptide", 5, min = 1, max=100, step = 1),
                       
                       # max covariance accepted
                       numericInput("MedCV_thresh", "Input maximal Covariance accepted", 0.65 , min = 0, max=1, step = 0.01),
                       
                       # remove peptides with missing data
                       numericInput("pNA", "Input missingness percentage threshold for peptide data", 0.99 , min = 0, max=1, step = 0.01)
      ),
      # choose log transform base
      selectInput("transform_base", "Choose method for protein data transformation", choices = c("log2", "log10", "sqrt", "quadratic", "BoxCox", "None")),
      # normalization method 
      selectInput("norm_method", "Choose method for protein data normalization", choices = c("SCoPE2", "None", "CONSTANd")),
      # multiple batch handling
      switchInput(inputId = "opts_multiple_batches", onLabel = "Advanced", offLabel = "Default", value = F, label="Options for multiple batches"),
      conditionalPanel(condition = "input.opts_multiple_batches", 
                       selectInput(inputId = "missing_v", "Choose method for missing value handling", choices=c("KNN", "drop rows", "replace with mean", "replace with median")),
                       selectInput(inputId = "batch_c", "Choose method for batch correction", choices=c("ComBat", "None"))
                       )
      
    ),
    
    # Main Panel for the output
    mainPanel(
      #create tabs
      tabsetPanel(type = "tabs",
                  tabPanel("Overview", plotOutput("overview_plot")),
                  tabPanel("Summary Barplot", plotOutput("summary_bar")),
                  tabPanel("Reporter Ion Intensity", 
                           selectInput("color_variable_ri", "select variable to indicate", ""),
                           plotOutput("RI_intensity")),
                  tabPanel("Covariance across razor peptides",
                           selectInput("color_variable_cv", "select variable to indicate", ""),
                           plotOutput("CV_median")),
                  tabPanel("Dimensionality reduction",
                           fluidRow(width=12,
                                    column(width=4,
                                           selectInput("color_variable_dim_red", "select variable to color", "")),
                                    column(width = 4,
                                           selectInput("shape_variable_dim_red", "select variable to shape", "")),
                                    column(width=4,
                                           selectInput("size_variable_dim_red", "select variable to size", ""))
                                    ),
                           fluidPage(
                           plotOutput("PCA"), 
                           plotOutput("UMAP"))
                           ),
                  tabPanel("Feature wise output", 
                           selectInput("selectedProtein", "Choose protein for observation", ""),
                           plotOutput("feature_subset")),
                  
                  tabPanel("Statistics",
                           actionButton("run_statistics", "Press to run statistics"),
                           actionButton("qqplot", "check dependencies for linear model"),
                           selectInput("model_design", "Choose your study design",
                                       choices = c("All pairwise comparison",
                                                   "Differential Expression with defined Contrasts",
                                                   "Multi factor additivity")),
                           uiOutput("sample_select"),
                           uiOutput("add_factor"),
                           plotOutput("venn_diagram"),
                           fluidRow(width=12,
                                    column(width=4,
                                           selectInput("chosen_coef", "Select your coefficient of interest", "")),
                                    column(width = 4,
                                           numericInput("p_value_cutoff", "choose p-value cutoff", value = 0.05, min = 0, step = 0.01)),
                                    column(width=4,
                                           numericInput("fold_change_cutoff", "choose fold change cutoff", value = 0.05, min = 0, step = 0.01))
                           ),
                           
                           plotlyOutput("volcano"),
                           tableOutput("protein_table")
                  )
      )
    )
  )
)

# box cox transformati0n --> alteration of basic MASS package was needed to run it
boxcox_1 <- function(object, ...) UseMethod("boxcox_1")

boxcox_1.formula <-
  function(object, lambda = seq(-2, 2, 1/10), eps=1/50)
  {
    object$y <- object$model
    m <- length(lambda)
    object <- stats::lm(object, y = TRUE, qr = TRUE)
    result <- NextMethod()
    result
  }

boxcox_1.lm <-
  function(object, lambda = seq(-2, 2, 1/10), eps=1/50)
  {
    object$y <- object$model
    m <- length(lambda)
    if(is.null(object$y) || is.null(object$qr))
      object <- update(object, y = TRUE, qr = TRUE)
    result <- NextMethod()
    result
  }

boxcox_1.default <-
  function(object, lambda = seq(-2, 2, 1/10), eps=1/50)
  {
    object$y <- as.matrix(object$model)
    if(is.null(y <- object$y) || is.null(xqr <- object$qr))
      stop(gettextf("%s does not have both 'qr' and 'y' components",
                    sQuote(deparse(substitute(object)))), domain = NA)
    if(any(y <= 0))
      stop("response variable must be positive")
    n <- length(y)
    ## scale y[]  {for accuracy in  y^la - 1 }:
    y <- y / exp(mean(log(y)))
    logy <- log(y) # now  ydot = exp(mean(log(y))) == 1
    xl <- loglik <- as.vector(lambda)
    m <- length(xl)
    for(i in 1L:m) {
      if(abs(la <- xl[i]) > eps)
        yt <- (y^la - 1)/la
      else yt <- logy * (1 + (la * logy)/2 *
                           (1 + (la * logy)/3 * (1 + (la * logy)/4)))
      loglik[i] <- - n/2 * log(sum(qr.resid(xqr, yt)^2))
    }
    list(x = xl, y = loglik)
  }


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)

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
  
  ## Analysis Pipeline ##
  # here is the basic processing from PSM to Protein done  
  
  # Data structures and objects
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
      
      evidence_data <- read.delim(input$evidence_file$datapath)
      meta_data_0 <- meta_data()

      if (input$file_level == FALSE) {
        
        incProgress(1/17, detail = paste("read Data"))
        scp_0 <- readSCP(featureData = evidence_data,
                         colData = meta_data_0,
                         channelCol = "Channel",
                         batchCol = "Raw.file",
                         removeEmptyCols = TRUE)
        
        #manual selection of SampleTypes to exclude from the analysis
        scp_0 <- scp_0[, !(scp_0$SampleType %in%  input$selectedSampleType_to_exclude)]      
        
        
        
        # # change zeros to NA, apply first filter
        incProgress(2/17, detail = paste("replacing zeros with NA"))
        scp_0 <- zeroIsNA(scp_0, 1:length(rowDataNames(scp_0)))
        
        # filter PSM
        # filter out potential contaminants
        # filter out matches to decoy database
        # keep PSMs with high PIF (parental ion fraction)
        incProgress(3/17, detail = paste("filter contaminants and PIF"))
        req(input$PIF_cutoff)
        PIF_cutoff<-input$PIF_cutoff
        scp_0 <- filterFeatures(scp_0,
                                ~ Reverse != "+" &
                                  Potential.contaminant != "+" &
                                  !is.na(PIF) & PIF > PIF_cutoff)
        
        # compute qvalues_PSMs to filter out by FDR
        incProgress(4/17, detail=paste("calculate q-value for PSMs"))
        scp_0 <- pep2qvalue(scp_0,
                            i = names(scp_0),
                            PEP = "PEP", # by reference the dart_PEP value is used
                            rowDataName = "qvalue_PSMs")
        
        incProgress(5/17, detail=paste("calculate q value for proteins"))
        scp_0 <- pep2qvalue(scp_0,
                            i = names(scp_0),
                            PEP = "PEP",
                            groupBy = "Leading.razor.protein",
                            rowDataName = "qvalue_proteins")
        
        incProgress(6/17, detail=paste("filter according to q-value"))      
        req(input$qvalue_cutoff)
        qvalue_cutoff <- input$qvalue_cutoff
        scp_0 <- filterFeatures(scp_0, ~ qvalue_proteins < qvalue_cutoff)
        
        # aggregate PSMS to peptides
        incProgress(7/17, detail=paste("aggregating features"))
        scp_0 <- aggregateFeaturesOverAssays(scp_0,
                                             i = names(scp_0),
                                             fcol = "Modified.sequence",
                                             name = paste0("peptides_", names(scp_0)),
                                             fun = matrixStats::colMedians, na.rm = TRUE)
        
        # join assays to one
        if (length(unique(meta_data_0$Raw.file)) > 1) {
          scp_0 <- joinAssays(scp_0,
                              i = ((length(names(scp_0))/2)+1):length(names(scp_0)),
                              name = "peptides")
        }
        
        
        # calculate median reporter IO intensity
        incProgress(8/17, detail=paste("calculate reporter ion intensity"))
        file_name <- unique(meta_data_0$Raw.file)
        peptide_file <- paste("peptides_", as.character(file_name), sep = "")
        
        if (length(peptide_file) > 1) {
          medians <- c()
          for (assay_name in peptide_file) {
            new_medians <- colMedians(assay(scp_0[[assay_name]]), na.rm = TRUE)
            medians <- c(medians, new_medians)
          }
        } else {
          medians <- colMedians(assay(scp_0[[peptide_file]]), na.rm = TRUE)
        }
        scp_0$MedianRI <- medians
        
        
        # Filter based on the median CV -> remove covariant peptides over multiple proteins
        incProgress(9/17, detail=paste("calculate covariance per cell"))
        req(input$nObs_pep_razrpr)
        nObs_pep_razrpr<-input$nObs_pep_razrpr
        
        if (length(peptide_file) > 1) {
          scp_0 <- medianCVperCell(scp_0,
                                   i = "peptides",
                                   groupBy = "Leading.razor.protein",
                                   nobs = nObs_pep_razrpr, 
                                   norm = "div.median",
                                   na.rm = TRUE,
                                   colDataName = "MedianCV")
        } else {
          scp_0 <- medianCVperCell(scp_0,
                                   i = peptide_file,
                                   groupBy = "Leading.razor.protein",
                                   nobs = nObs_pep_razrpr, 
                                   norm = "div.median",
                                   na.rm = TRUE,
                                   colDataName = "MedianCV")
        }
        
        incProgress(10/17, detail=paste("filtering according to covariance"))
        req(input$MedCV_thresh)
        MedCV_thresh <- input$MedCV_thresh
        scp_0 <- scp_0[, !is.na(scp_0$MedianCV) & scp_0$MedianCV < MedCV_thresh, ]
        
        incProgress(12/17, detail=paste("remove peptides by missing rate"))
        req(input$pNA)
        pNA <- input$pNA
        
        if (length(peptide_file) > 1) {
          scp_0 <- filterNA(scp_0,
                            i = "peptides",
                            pNA = pNA)
        } else {
          scp_0 <- filterNA(scp_0,
                            i = peptide_file,
                            pNA = pNA)
        }
        
        incProgress(13/17, detail=paste("aggregate peptides to proteins"))
        
        # aggregate peptide to protein
        if (length(peptide_file) > 1) {
          scp_0 <- aggregateFeatures(scp_0,
                                     i = "peptides",
                                     name = "proteins",
                                     fcol = "Leading.razor.protein",
                                     fun = matrixStats::colMedians, na.rm = TRUE)
          
        } else {
          scp_0 <- aggregateFeatures(scp_0,
                                     i = peptide_file,
                                     name = "proteins",
                                     fcol = "Leading.razor.protein",
                                     fun = matrixStats::colMedians, na.rm = TRUE)
        }        
        
        
      } 

      
      incProgress(14/17, detail=paste("transforming protein data"))
      req(input$transform_base)

      if (input$transform_base == "log2") {
        if (input$file_level == FALSE) {
          scp_0 <- logTransform(scp_0,
                                base = 2,
                                i = "proteins",
                                name = "proteins_transf")
        } else {
          evidence_data <- log2(evidence_data)
        }
      }
      else if (input$transform_base == "log10") {
        if (input$file_level == FALSE) {
          scp_0 <- logTransform(scp_0,
                                base = 10,
                                i = "proteins",
                                name = "proteins_transf")
        } else {
          evidence_data <- log10(evidence_data)
        }
      }
      else if (input$transform_base == "sqrt") {
        if (input$file_level == FALSE) {
          scp_0 <- sweep(scp_0, i="proteins",
                         MARGIN = 2,
                         FUN="^",
                         STATS=1/2,
                         name="proteins_transf")
        } else {
          evidence_data <- sqrt(evidence_data)
        }
      }
      else if (input$transform_base == "quadratic") {
        if (input$file_level == FALSE) {
          scp_0 <- sweep(scp_0, i="proteins",
                         MARGIN = 2,
                         FUN="^",
                         STATS=2,
                         name="proteins_transf")
        } else {
          evidence_data <- evidence_data**2
        }
      }
      else if (input$transform_base == "BoxCox") {
        if (input$file_level == FALSE) {
          protein_matrix <- assay(scp_0[["proteins"]])
          b <- boxcox_1(stats::lm(protein_matrix ~ 1))
          # Exact lambda
          lambda <- b$x[which.max(b$y)]
          print(lambda)
          
          
          if (round(lambda, digits = 0) == -2 || lambda < -1.5) {
            protein_matrix <- 1/protein_matrix**2
            print("1/protein_matrix**2")
          }
          if (round(lambda, digits = 0) == -1 || lambda < -0.75 && lambda > -1.5) {
            protein_matrix <- 1/protein_matrix
            print("1/protein_matrix")
          }
          if (round(lambda, digits = 1) == -0.5 || lambda < -0.25 && lambda > -0.75) {
            protein_matrix <- 1/(protein_matrix**1/2)
            print("1/(protein_matrix**1/2)")
          }
          if (round(lambda, digits = 0) == 0 || lambda < 0.25 && lambda > - 0.25 ) {
            protein_matrix <- log10(protein_matrix)
            transform_base_bc <- "log10"
            print("log10(protein_matrix)")
          }
          if (round(lambda, digits = 1) == 0.5 || lambda > 0.25 && lambda < 0.75) {
            protein_matrix <- protein_matrix**1/2
            print("protein_matrix**1/2")
          }
          if (round(lambda, digits = 0) == 1 || lambda > 0.75 && lamdba < 1.5) {
            protein_matrix <- protein_matrix
            print("protein_matrix")
          }
          if (round(lambda, digits = 0) == 2 || lambda > 1.5) {
            protein_matrix <- protein_matrix**2
            print("protein_matrix**2")
          }        
          
          sce <- getWithColData(scp_0, "proteins")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_transf")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins",
                                        to = "proteins_transf")
          
          assay(scp_0[["proteins_transf"]]) <- protein_matrix
        } else {
          protein_matrix <- evidence_data
          b <- boxcox_1(stats::lm(evidence_data ~ 1))
          # Exact lambda
          lambda <- b$x[which.max(b$y)]
          print(lambda)
          
          if (round(lambda, digits = 0) == -2 || lambda < -1.5) {
            protein_matrix <- 1/protein_matrix**2
            print("1/protein_matrix**2")
          }
          if (round(lambda, digits = 0) == -1 || lambda < -0.75 && lambda > -1.5) {
            protein_matrix <- 1/protein_matrix
            print("1/protein_matrix")
          }
          if (round(lambda, digits = 1) == -0.5 || lambda < -0.25 && lambda > -0.75) {
            protein_matrix <- 1/(protein_matrix**1/2)
            print("1/(protein_matrix**1/2)")
          }
          if (round(lambda, digits = 0) == 0 || lambda < 0.25 && lambda > - 0.25 ) {
            protein_matrix <- log10(protein_matrix)
            transform_base_bc <- "log10"
            print("log10(protein_matrix)")
          }
          if (round(lambda, digits = 1) == 0.5 || lambda > 0.25 && lambda < 0.75) {
            protein_matrix <- protein_matrix**1/2
            print("protein_matrix**1/2")
          }
          if (round(lambda, digits = 0) == 1 || lambda > 0.75 && lamdba < 1.5) {
            protein_matrix <- protein_matrix
            print("protein_matrix")
          }
          if (round(lambda, digits = 0) == 2 || lambda > 1.5) {
            protein_matrix <- protein_matrix**2
            print("protein_matrix**2")
          }        
          evidence_data <- protein_matrix
        }
      }
      else if (input$transform_base == "None") {
        if (input$file_level == FALSE) {
          sce <- getWithColData(scp_0, "proteins")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_transf")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins",
                                        to = "proteins_transf")        
        } else {
          evidence_data <- evidence_data
        }
      }
      
      incProgress(15/17, detail=paste("normalizing proteins"))
      req(input$norm_method)
      if (input$norm_method == "SCoPE2" && input$transform_base == "log2" | input$transform_base == "log10") {
        if (input$file_level == FALSE) {
          # center cols with median
          scp_0 <- sweep(scp_0, i = "proteins_transf",
                         MARGIN = 2,
                         FUN = "-",
                         STATS = colMedians(assay(scp_0[["proteins_transf"]]),
                                            na.rm = TRUE),
                         name = "proteins_norm_col")
          
          # Center rows with mean
          scp_0 <- sweep(scp_0, i = "proteins_norm_col",
                         MARGIN = 1,
                         FUN = "-",
                         STATS = rowMeans(assay(scp_0[["proteins_norm_col"]]),
                                          na.rm = TRUE),
                         name = "proteins_norm")
          
        } else {
          # #normalize colwise
          evidence_data <- sweep(evidence_data[,-1], 2, colMedians(evidence_data[,-1]), FUN="-")
          # #normalize rowwise
          evidence_data <- sweep(evidence_data[,-1], 1, rowMeans(evidence_data[,-1]), FUN="-")
        }
        
      } else if (input$norm_method == "SCoPE2" && input$transform_base != "log2" | input$transform_base != "log10") {
        if (input$file_level == FALSE) {
          # center cols with median
          scp_0 <- sweep(scp_0, i = "proteins_transf",
                         MARGIN = 2,
                         FUN = "/",
                         STATS = colMedians(assay(scp_0[["proteins_transf"]]),
                                            na.rm = TRUE),
                         name = "proteins_norm_col")
          
          # Center rows with mean
          scp_0 <- sweep(scp_0, i = "proteins_norm_col",
                         MARGIN = 1,
                         FUN = "/",
                         STATS = rowMeans(assay(scp_0[["proteins_norm_col"]]),
                                          na.rm = TRUE),
                         name = "proteins_norm")
          
        } else {
          # #normalize colwise
          evidence_data <- sweep(evidence_data[,-1], 2, colMedians(evidence_data[,-1]), FUN="/")
          # #normalize rowwise
          evidence_data <- sweep(evidence_data[,-1], 1, rowMeans(evidence_data[,-1]), FUN="/")
        }
        
      } else if (input$norm_method == "CONSTANd") {
        if (input$file_level == FALSE) {
          # apply matrix raking --> row means and col means equal Nrows and Ncols  
          protein_matrix <- assay(scp_0[["proteins_transf"]])
          protein_matrix <- CONSTANd(protein_matrix)
          
          sce <- getWithColData(scp_0, "proteins_transf")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_norm")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins_transf",
                                        to = "proteins_norm")
          
          assay(scp_0[["proteins_norm"]]) <- protein_matrix$normalized_data
        } else {
          evidence_data <- CONSTANd(evidence_data)
        }
      } else if (input$norm_method == "None") {
        if (input$file_level == FALSE) {
          sce <- getWithColData(scp_0, "proteins_transf")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_norm")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins_transf",
                                        to = "proteins_norm")
        } else {
          evidence_data <- evidence_data
        }
      }
      
      
      if (length(peptide_file) > 1) {
        incProgress(16/17, detail=paste("running missing value imputation"))
        if (input$missing_v == "KNN") {
          scp_0 <- impute(scp_0,
                          i = "proteins_norm",
                          name = "proteins_imptd",
                          method = "knn",
                          k = 3, rowmax = 1, colmax= 1,
                          maxp = Inf, rng.seed = as.numeric(gsub('[^0-9]', '', Sys.Date())))
        } else if (input$missing_v == "drop rows") {
          sce <- getWithColData(scp_0, "proteins_norm")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_imptd")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins_norm",
                                        to = "proteins_imptd")
          
          scp_0 <- filterNA(scp_0, pNA = 0, "proteins_imptd")
        } else if (input$missing_v == "replace with mean") {
          sce <- getWithColData(scp_0, "proteins_norm")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_imptd")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins_norm",
                                        to = "proteins_imptd")
          
          assay(scp_0[["proteins_imptd"]]) <- replace(assay(scp_0[["proteins_imptd"]]), is.na(assay(scp_0[["proteins_imptd"]])), mean(assay(scp_0[["proteins_imptd"]]), na.rm = TRUE))
        } else if (input$missing_v == "replace with median") {
          sce <- getWithColData(scp_0, "proteins_norm")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_imptd")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins_norm",
                                        to = "proteins_imptd")
          
          assay(scp_0[["proteins_imptd"]]) <- replace(assay(scp_0[["proteins_imptd"]]), is.na(assay(scp_0[["proteins_imptd"]])), median(assay(scp_0[["proteins_imptd"]]), na.rm = TRUE))        
        }
        
        incProgress(16/17, detail=paste("running batch correction"))
        sce <- getWithColData(scp_0, "proteins_imptd")
        if (input$batch_c == "ComBat") {
          batch <- colData(sce)$Raw.file
          # can be used to aim batch correction for desired result
          model <- model.matrix(~0 + SampleType, data = colData(sce))
          
          assay(sce) <- ComBat(dat = assay(sce),
                               batch = batch)#,
          #mod = model)
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_batchC")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins_imptd",
                                        to = "proteins_batchC")
          
          sce <- getWithColData(scp_0, "proteins_batchC")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_dim_red")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins_batchC",
                                        to = "proteins_dim_red") 
          
        } else if (input$batch_c == "None") {
          sce <- getWithColData(scp_0, "proteins_imptd")
          
          scp_0 <- addAssay(scp_0,
                            y = sce,
                            name = "proteins_dim_red")
          
          scp_0 <- addAssayLinkOneToOne(scp_0,
                                        from = "proteins_imptd",
                                        to = "proteins_dim_red")          
        }
        
        
        
        
      } else {
        sce <- getWithColData(scp_0, "proteins_norm")
        
        scp_0 <- addAssay(scp_0,
                          y = sce,
                          name = "proteins_dim_red")
        
        scp_0 <- addAssayLinkOneToOne(scp_0,
                                      from = "proteins_norm",
                                      to = "proteins_dim_red")  
      }
      
      
      
      
      incProgress(16/17, detail=paste("running dimensionality reduction"))
      
      scp_0[["proteins_dim_red"]] <- scater::runPCA(scp_0[["proteins_dim_red"]],
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
  
  # expression matrix is used by plot functions
  exp_matrix <- reactive({
    scp_0 <- scp()
    exp_matrix_0 <- data.frame(assay(scp_0[["proteins_dim_red"]]))
    colnames(exp_matrix_0) <- scp_0$SampleType  
    return(exp_matrix_0)
  })
  
  # trigger for the CONSTANd dependency MA plot
  CONSTANd_trigger <- reactive(list(input$update_button, input$norm_method))
  ### graphical outputs and user interface
  
  
  # observer for protein list
  observe({
    req(protein_list())
    updateSelectInput(session, "selectedProtein", choices = protein_list())
  })
  
  # reactive element for the protein list --> update if the scp object changes
  protein_list <- reactive({
    req(scp())
    req(scp()[["proteins"]])
    scp_0 <- scp()
    list <- rowData(scp_0)[["proteins"]][,1]
    return(list)
  })
  

  # observer for the CONSTANd normalization dependency to check before analysis
  observeEvent(ignoreInit=TRUE, CONSTANd_trigger(), {
    if (CONSTANd_trigger()[2] == "CONSTANd" && !is.null(comp_list())){
      validate(need(input$update_button > 0, ''))
      showModal(CONSTANdModal())
    }
  })
  
  # observer for comparison between sample types
  observeEvent(CONSTANd_trigger(), {
    updateSelectInput(session, "selectedComp", choices = comp_list())
  })
  
  
  
  # reactive element for the comp list --> update if the scp object changes
  comp_list <- reactive({
    req(scp())
    scp_0 <- scp()
    list <- apply(combn(unique(scp_0$SampleType),2),2,paste, collapse="-")
    return(list)
  })

  columns <- reactive({
    req(scp())
    scp_0 <- scp()
    list <- colnames(colData(scp_0))
    return(list)
  })
  
  # observer for sample type exclusion
  observe({
    req(sample_types())
    updateSelectInput(session, "selectedSampleType_to_exclude", choices = sample_types())
  })
  sample_types <- reactive({
    req(scp())
    scp_0 <- scp()
    list <- scp_0$SampleType
    return(list)
  })
  
  ## Output of analysis pipeline 
  # pathway of the data first tab
  output$overview_plot <- renderPlot({
    if (!is.null(scp())) {
      plot(scp()) }
  })
  
  # summary barchart second tab
  output$summary_bar <- renderPlot({
    scp_0 <- scp()

    count_table <- data.frame(dims(scp_0))
    
    nPSMs <- log10(sum(count_table[1:length(unique(scp_0$Raw.file))]))
    nPeps <- log10(sum(count_table$peptides))
    nProts <- log10(sum(count_table$proteins))
    summary_count <- data.frame(nPSMs, nPeps, nProts)
    summary_count <- reshape2::melt(summary_count)
    
    ggplot(summary_count, aes(x = variable, y = value, fill = variable)) +
      geom_bar(stat = "identity") +
      labs(x = "log10 of number of PSMs, Peptides and Proteins", y = "Counts") +
      scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"))
  })
  
  # Reporter Ion intensity visualisation
  #observer for color_variable
  observe({
    req(columns())
    updateSelectInput(session, "color_variable_ri", choices=columns())
  })  
  # Boxplot of reporter ion intensity third tab
  output$RI_intensity <- renderPlot({
    scp_0 <- scp()
    
    colData(scp_0) %>%
      data.frame %>%
      ggplot() +
      aes(x = MedianRI, 
          y = get(input$color_variable_ri),
          fill = get(input$color_variable_ri)) +
      geom_boxplot() +
      scale_x_log10() +
      labs(fill=as.character(input$color_variable_ri), y=as.character(input$color_variable_ri)) 
  })
  
  # Covariance visualisation
  #observer for color_variable
  observe({
    req(columns())
    updateSelectInput(session, "color_variable_cv", choices=columns())
  })  
  # covariance across razor peptides fourth tab
  output$CV_median <- renderPlot({
    req(scp())
    req(meta_data())
    
    scp_0 <- scp()
    meta_data_0 <- meta_data()
    
    file_name <- unique(meta_data_0$Raw.file)
    peptide_file <- paste("peptides_", as.character(file_name), sep = "")
    
    if (length(peptide_file) > 1) {
      getWithColData(scp_0, "peptides") %>%
        colData %>%
        data.frame %>%
        ggplot(aes(x = MedianCV,
                   fill = get(input$color_variable_cv))) +
        geom_boxplot()+
        labs(fill=as.character(input$color_variable_cv), y=as.character(input$color_variable_cv))  
    } else {
      getWithColData(scp_0, peptide_file) %>%
        colData %>%
        data.frame %>%
        ggplot(aes(x = MedianCV,
                   fill = get(input$color_variable_cv))) +
        geom_boxplot()+
        labs(fill=as.character(input$color_variable_cv), y=as.character(input$color_variable_cv))
    }
  })
  
  ## Dimensionality reduction plot
  #observer for color_variable
  observe({
    req(columns())
    updateSelectInput(session, "color_variable_dim_red", choices=c("NULL", columns()))
  })
  #observer for shape_variable
  observe({
    req(columns())
    updateSelectInput(session, "shape_variable_dim_red", choices=c("NULL", columns()))
  })
  #observer for size_variable
  observe({
    req(columns())
    updateSelectInput(session, "size_variable_dim_red", choices=c("NULL", columns()))
  })
  
  # principle component analysis in fifth tab
  output$PCA <- renderPlot({
    extract_null <- function(variable) if (variable == "NULL") {
      return(NULL)
    } else {
      return(variable)
    }
    
    color_variable_dim_red <- extract_null(input$color_variable_dim_red)
    shape_variable_dim_red <- extract_null(input$shape_variable_dim_red)
    size_variable_dim_red <- extract_null(input$size_variable_dim_red)
    
    scp_0 <- scp()
    plotReducedDim(scp_0[["proteins_dim_red"]],
                   dimred = "PCA",
                   colour_by = color_variable_dim_red,
                   shape_by = shape_variable_dim_red,
                   size_by = size_variable_dim_red,
                   point_alpha = 1,
                   point_size=3)
  })
  
  # Umap dimensionality reduction in fith tab
  output$UMAP <- renderPlot({
    extract_null <- function(variable) if (variable == "NULL") {
      return(NULL)
    } else {
      return(variable)
    }
    
    color_variable_dim_red <- extract_null(input$color_variable_dim_red)
    shape_variable_dim_red <- extract_null(input$shape_variable_dim_red)
    size_variable_dim_red <- extract_null(input$size_variable_dim_red)
    
    
    
    scp_0 <- scp()
    plotReducedDim(scp_0[["proteins_dim_red"]],
                   dimred = "UMAP",
                   colour_by = color_variable_dim_red,
                   shape_by = shape_variable_dim_red,
                   size_by = size_variable_dim_red,
                   point_alpha = 1,
                   point_size = 3)
  })
  
  
  # protein wise visualisation 
  output$feature_subset <- renderPlot({
    scp_0 <- scp()
    subsetByFeature(scp_0, input$selectedProtein) %>%
      ## Format the `QFeatures` to a long format table
      longFormat(colvars = c("Raw.file", "SampleType", "Channel")) %>%
      data.frame %>%
      ## This is used to preserve ordering of the samples and assays in ggplot2
      mutate(assay = factor(assay, levels = names(scp_0)),
             Channel = sub("Reporter.intensity.", "", Channel)) %>%
      mutate(Channel = as.numeric(Channel)) %>%
      arrange(Channel) %>%
      mutate(Channel = factor(Channel, levels = unique(Channel))) %>%
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

  # plot MA plot for the constand method
  output$MAplot <- renderPlot({
    req(input$selectedComp)
    req(scp())
    
    scp_0 <- scp() 
    scp_0 <- filterNA(scp_0, pNA = 0, "proteins_transf")
    
    user_choice <- input$selectedComp
    # split user choice of comp back to sample types
    user_choice_vector <- strsplit(user_choice, split = "-")
    # and assign them to a variable
    choice_A <- user_choice_vector[[1]][1]
    choice_B <- user_choice_vector[[1]][2]
    
    #find the row indeces of corresponding to the individual sample types
    st_indeces <- split(seq_along(scp_0$SampleType), scp_0$SampleType)
    index_A <- st_indeces[choice_A]
    index_B <- st_indeces[choice_B]
    
    
    MAplot <- function(x,y,use.order=FALSE, R=NULL, cex=1.6, showavg=TRUE) {
      # catch unequal size of matrices
      if (dim(x)[2] != dim(y)[2]) {
        if (dim(x)[2] > dim(y)[2]) {
          x <- x[, 1:dim(y)[2]]
        } else if (dim(y)[2] > dim(x)[2]) {
          y <- y[, 1:dim(x)[2]]
        }
      }
      
      # make an MA plot of y vs. x that shows the rolling average,
      M <- log2(y/x)
      xlab = 'A'
      if (!is.null(R)) {r <- R; xlab = "A (re-scaled)"} else r <- 1
      A <- (log2(y/r)+log2(x/r))/2
      if (use.order) {
        orig.order <- order(A)
        A <- orig.order
        M <- M[orig.order]
        xlab = "original rank of feature magnitude within IPS"
      }
      # select only finite values
      use <- is.finite(M)
      A <- A[use]
      M <- M[use]
      # plot
      print(var(M))
      plot(A, M, xlab=xlab, cex.lab=cex, cex.axis=cex)
      # rolling average
      if (showavg) { lines(lowess(M~A), col='red', lwd=5) }
    }
    MAplot(assay(scp_0[["proteins_transf"]][,index_A[[1]]]), assay(scp_0[["proteins_transf"]][,index_B[[1]]]))
  })
  
  #create interface for the constand normalization method
  CONSTANdModal <- function() {
    modalDialog(
      selectInput("selectedComp", "Choose comparison for observation", ""),
      plotOutput("MAplot"),
      HTML(paste("
        <h2> Data assumptions </h2>
        <p>To warrant that the normalization procedure itself is not biased, three assumptions about the data are made (as for any type of data-driven normalization method) which may be verified by making MA-plots.
        <ul>
        <li><strong>The majority of features (peptides/genes/…) are not differentially expressed</strong>, to avoid a bias in the estimate of the mean value. It is assumed that up-down shifts are not due to biological causes. The reference set used in the normalization step is the set of all peptides identified in the experiment.<br>
        <em>MA-plot: the observations form a single ‘cloud’ with a dense center and less dense edges, as opposed to, for instance, two clouds or a cloud with uniform density</em>.</li>
        <li><strong>The number of up-regulated features is roughly equal to the number of down-regulated ones</strong>. If the data were skewed, this would lead to a bias in the normalization result.<br>
        <em>MA-plot: the cloud of observations exhibits a bilateral symmetry about some axis (usually horizontal, but inclinations may occur)</em>.</li>
        <li><strong>Any systematic bias present is linearly proportional to the magnitude of the quantification values</strong>. Only this way it is possible to find one appropriate normalization factor for each quantification sample.<br>
        <em>MA-plot: the axis of bilateral symmetry is a straight line (which may be inclined), i.e., the moving average M-values of the central cloud form an approximately <strong>straight</strong> line</em>.</li>
        </ul>
      ")),
      footer = tagList(
        modalButton("Dismiss")
      )
    )
  }
  
  
    

  ### Statistic Module ###
  observeEvent(input$model_design, {
    updateSelectInput(session, "selectedComp_stat", choices = sample_types())
  })
  
  observeEvent(input$model_design, {
    updateSelectInput(session, "col_factors", choices = columns())
  })
  
  #observer for updating the coefficient dropdown menu
  observeEvent(input$run_statistics, {
    updateSelectInput(session, "chosen_coef", choices = coefs())
  })
  
  
  # statistical pipeline
  stat_result <- eventReactive(input$run_statistics, {
    withProgress(message= "running statistical analysis", value=0, {
      incProgress(1/3, detail=paste("read data"))
      scp_0 <- scp()
      exp_matrix_0 <- exp_matrix()
      
      incProgress(2/3, detail=paste("creating linear model"))
      req(input$model_design)
      # Create a design matrix
      if (input$model_design == "All pairwise comparison") {
        design <- model.matrix(~0+factor(scp_0$SampleType))
        colnames(design) <- unique(scp_0$SampleType)
        fit <- lmFit(exp_matrix_0, design)
      }
      #Differential Expression with defined Contrasts
      else if (input$model_design == "Differential Expression with defined Contrasts") {
        # fetch user selection
        req(input$selectedComp_stat)
        scp_0 <- scp_0[, scp_0$SampleType %in%  input$selectedComp_stat]
        exp_matrix_0 <- data.frame(assay(scp_0[["proteins_dim_red"]]))
        colnames(exp_matrix_0) <- scp_0$SampleType
        
        design <- model.matrix(~0+factor(scp_0$SampleType))
        colnames(design) <- unique(scp_0$SampleType)
        
        fit <- lmFit(exp_matrix_0, design)
        
        user_contrast <- paste(input$selectedComp_stat, sep = "-", collapse = NULL)
        cont_matrix <- makeContrasts(contrasts=user_contrast,levels=colnames(design))
        fit <- contrasts.fit(fit, cont_matrix)
      }
      else if (input$model_design == "Multi factor additivity") {
        req(input$col_factors)
        req(input$selectedComp_stat)
        
        scp_0 <- scp_0[, scp_0$SampleType %in%  input$selectedComp_stat]
        exp_matrix_0 <- data.frame(assay(scp_0[["proteins_dim_red"]]))
        colnames(exp_matrix_0) <- scp_0$SampleType
        
        fetched_factor <- colData(scp_0)[input$col_factors]
        design_frame <- cbind(fetched_factor, scp_0$SampleType)
        design <- model.matrix(~0+ . , data=design_frame)
        fit <- lmFit(exp_matrix_0, design)
      }
      
      incProgress(4/4, detail=paste("Bayes statistics of differential expression"))
      # *There are several options to tweak!*
      fit <- eBayes(fit)
    })
    return(fit)
  })

  ## Dependencies  
  # initial qq counter for the first plot 
  qq_count <- reactiveVal(1)
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
  # plot the qq_count-th column
  output$qqPlot <- renderPlot({
    exp_matrix_0 <- exp_matrix()
    sample_types <- unique(colnames(exp_matrix_0))
    i <- qq_count()
    qqnorm(exp_matrix_0[[sample_types[i]]], main=paste("QQ-plot of ", sample_types[i]))
    qqline(exp_matrix_0[[sample_types[i]]])
  })
  # plot hist for linear model
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

  # additional ui for statistic module
  output$sample_select <- renderUI({
    req(input$model_design == "Differential Expression with defined Contrasts" | input$model_design == "Multi factor additivity")
    selectInput("selectedComp_stat", "choose your samples of interest", "", multiple = T)
  })

  # Ui for multifactorial model
  output$add_factor <- renderUI({
    req(input$model_design == "Multi factor additivity")
    selectInput("col_factors", "choose second factor or multiple for your model", "", multiple = T)
  })
  
  
  ## Results 
  # coefficients for the statistic modules volcanoplot and topTable
  coefs <- reactive({
    req(stat_result())
    fit_0 <- stat_result()
    list <- colnames(fit_0$coefficients)
  })
  
  # volcanoplot in statistic tab
  output$volcano <- renderPlotly({
    # load data.frame
    toptable <- protein_table()
    toptable <- cbind(rownames(toptable), data.frame(toptable, row.names=NULL))
    colnames(toptable)[1] <- "protein"
    
    # significance cutoffs
    FCcutoff <- input$fold_change_cutoff
    pCutoff <- input$p_value_cutoff
    
    # define factor levels for significance
    toptable$Significance <- "not significant"
    toptable$Significance[(abs(toptable$logFC) > FCcutoff)] <- "fold change"
    toptable$Significance[(toptable$P.Value < pCutoff)] <- "p-value"
    toptable$Significance[(toptable$P.Value < pCutoff) & (abs(toptable$logFC) > 
                                                            FCcutoff)] <- "fold change & p-value"
    toptable$Significance <- factor(toptable$Sig, levels = c("not significant", "fold change", 
                                                             "p-value", "fold change & p-value"))
    
    p <- ggplot(toptable, aes(x = logFC, y = -log10(P.Value), text=protein)) +
      geom_point(aes(color = Significance), alpha = 1/2, shape = 19, size = 1.5, na.rm = TRUE) +
      theme(legend.position = "top")
    
    ggplotly(p, tooltip = "text")

  })
  
  # reactive element for toptable
  protein_table <- reactive({
    req(stat_result())
    req(input$chosen_coef)
    data.frame(topTable(stat_result(), number = Inf, adjust = "BH", coef = input$chosen_coef))
  })
  
  # show significant proteins in the table
  output$protein_table <- renderTable({
    req(protein_table())
    protein_table()
  }, rownames = T, striped = T)
  
  # venn diagram for significant proteins
  output$venn_diagram <- renderPlot({
    req(stat_result())
    vennDiagram(decideTests(stat_result()))
  })
}


# Run the application 
shinyApp(ui = ui, server = server)