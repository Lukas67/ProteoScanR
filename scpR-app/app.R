library("shiny")
library("shinyWidgets")
library("shinydashboard")
library("DT")
library("spsComps")

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
library("clusterProfiler")
library("infotheo")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("piano")
library("snowfall")
library("parallel")
library("visNetwork")



reactlog::reactlog_enable()

# Define UI for application that draws a histogram
ui <- fluidPage(
  # useShinyjs(),
  # Application title
  titlePanel(fluidRow(width=15,
                      column(h1("Proteomics Workbench"), h5("by Lukas Gamp"), offset = 0.5, width=5))
    ),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(id="settings_pane",
      # the button to rule them all
      actionButton("update_button", "Press to run/update"),
      actionButton("help", "Help"),
      
      # read in of the data
      materialSwitch(inputId = "file_level", value = F, label="Handle preprocessed expression set", status = "danger"),
      fileInput("evidence_file", "Upload tabdel (.txt) PSM file", accept = c("text")),
      fileInput("sample_annotation_file", "Upload tabdel (.txt) sample annotation file", accept = c("text")),
      selectInput("selectedSampleType_to_exclude", "Select SampleType to exclude", "", multiple=T),
      
      conditionalPanel(id="full_control_pane", condition = "!input.file_level", 
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
      selectInput("norm_method", "Choose method for protein data normalization", choices = c("col-Median,row_Mean", "None", "CONSTANd", "Quantile")),
      # multiple batch handling
      switchInput(inputId = "opts_multiple_batches", onLabel = "Advanced", offLabel = "Default", value = F, label="Options for multiple batches"),
      conditionalPanel(id="multiple_batch_pane",condition = "input.opts_multiple_batches", 
                       selectInput(inputId = "missing_v", "Choose method for missing value handling", choices=c("KNN", "drop rows", "replace with mean", "replace with median")),
                       selectInput(inputId = "batch_c", "Choose method for batch correction", choices=c("None", "ComBat"))
      )
      
    ),
    
    # Main Panel for the output
    mainPanel(id="output_pane",
      #create tabs
      tabsetPanel(id="outbut_tabset", type = "tabs",
                  tabPanel("Overview", id="overview_pane",
                           verbatimTextOutput("debug"),
                           plotOutput("overview_plot", width = "1500px", height = "1000px")),
                  tabPanel("Summary Barplot", id="summary_pane",
                           plotOutput("summary_bar", width = "1300px", height = "1000px")),
                  tabPanel("Normalization validation", id="normval_pane",
                           selectInput("selectedSampleType_normval", "Select SampleType to validate", ""),
                           plotOutput("norm_val_plot_intra_group"),
                           selectInput("selectedComp_normval", "Choose comparison for observation", ""),
                           plotOutput("norm_val_plot_inter_group")
                  ),
                  tabPanel("Reporter Ion Intensity", id="ri_pane",
                           selectInput("color_variable_ri", "select variable to indicate", ""),
                           plotOutput("RI_intensity", width = "1300px", height = "800px")),
                  tabPanel("Covariance and correlation", id="cov_cor_pane",
                           fluidRow(width=10,
                                    conditionalPanel(condition = "!input.file_level", id="razor_prot_pane",
                                                     column(width = 5, conditionalPanel(condition = "input.cv_plot_switch",
                                                                                        selectInput("color_variable_cv", "select variable to indicate", ""))
                                                                             ),

                                                     column(width = 5, switchInput("cv_plot_switch", label="switch between visualizations", onLabel="Cov boxlot", offLabel="Cor heatmap", onStatus="primary", offStatus = "info"))
                                                     )
                                    ),
                           conditionalPanel(condition = "!input.cv_plot_switch",
                                            plotOutput("corr_matrix", width = "1400px", height = "1000px")
                                            ),
                           conditionalPanel(condition = "input.cv_plot_switch",
                                            plotOutput("CV_median", width = "1200px", height = "900px")
                                            )
                  ),
                  tabPanel("Dimensionality reduction", id="dimred_pane",
                           materialSwitch(inputId = "third_dim", value = F, label="3d", status = "danger"),
                           fluidRow(width=12,
                                    column(width=4,
                                           selectInput("color_variable_dim_red", "select variable to color", "")),
                                    conditionalPanel(condition = "!input.third_dim", id="dimred_3d_set_pane",
                                                     column(width = 4,
                                                            selectInput("shape_variable_dim_red", "select variable to shape", "")),
                                                     column(width=4,
                                                            selectInput("size_variable_dim_red", "select variable to size", ""))
                                    )
                           ),
                           conditionalPanel(condition = "!input.third_dim", id="dimred_2d_pane",
                                            fluidPage(
                                              h4("PCA-plot"),
                                              plotOutput("PCA"),
                                              h4("UMAP-plot"),
                                              plotOutput("UMAP")
                                            )
                           ),
                           conditionalPanel(condition = "input.third_dim", id="dimred_3d_pane",
                                            fluidPage(
                                              h4("PCA-plot"),
                                              plotlyOutput("PCA_third_dim"),
                                              h4("UMAP-plot"),
                                              plotlyOutput("UMAP_third_dim")
                                            )
                           )
                  ),
                  tabPanel("Feature wise output", id="feature_wise_pane",
                           selectizeInput("selectedProtein", "Choose protein for observation", ""),
                           conditionalPanel(condition = "input.file_level", id="color_settings_pane",
                                            selectInput("color_variable_feature", "select variable to color", "")
                           ),
                           plotOutput("feature_subset", width = "1300px", height = "1000px")),
                  
                  tabPanel("Statistics", id="stats_pane",
                           actionButton("run_statistics", "Press to run statistics"),
                           actionButton("qqplot", "check dependencies for linear model"),
                           fluidRow(width=12,
                                    column(width=4,
                                           selectInput("model_design", "choose your study design",
                                                       choices = c("All pairwise comparison",
                                                                   "Differential Expression with defined Contrasts",
                                                                   "Multi factor additivity"))),
                                    column(width = 4,
                                           uiOutput("sample_select")),
                                    column(width=4,
                                           uiOutput("add_factor"))
                           ),
                           selectInput("p_value_correction", "select method for p-value correction", choices = c("none", "BH", "fdr", "BY", "holm")),
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
                           DT::dataTableOutput("protein_table")
                  ),
                  tabPanel("Protein set enrichment Pathway based", id="psa_pathw_pane",
                           br(),
                           actionButton("run_kegg", "Press to run pathway enrichment", class="btn-success"),
                           br(),
                           br(),
                           fluidRow(width=12,
                                    column(width = 3,
                                           switchInput("pw_plot_switch", label="switch between visualizations", onLabel="Bar/Dot", offLabel="Network", onStatus="primary", offStatus = "info")),
                                    column(width = 3,
                                           materialSwitch("design_plot_gsea", label = "change Plot design")),
                                    column(width = 3,
                                           materialSwitch("panel_as_universe", label = "use detected proteins as background")),
                                    column(width=3,
                                           fileInput("custom_universe", "Upload tabdel (.txt) uniprot_id background (w.col header)", accept = c("text")))
                           ),
                           br(),
                           fluidRow(width=12,
                                    column(width=4,
                                           textOutput("p_correct_1")),
                                    column(width = 4,
                                           textOutput("p_cutoff_1")),
                                    column(width=4,
                                           textOutput("fc_cutoff_1"))
                           ),
                           conditionalPanel("!input.pw_plot_switch", id="network_pane",
                                            conditionalPanel("input.design_plot_gsea",
                                                             visNetworkOutput("gsea_network_1", width = "1300px", height="800px")
                                                             ),
                                            conditionalPanel("!input.design_plot_gsea",
                                                             plotOutput("gsea_network_2", width = "1300px", height="800px")
                                                             )
                           ),
                           conditionalPanel("input.pw_plot_switch",
                                            plotlyOutput("gsea", width="1300px", height = "800px"),
                           ),
                
                  ),
                  tabPanel("Protein set enrichment Ontology based", id="psa_ont_pane",
                           br(),
                           actionButton("go_ontology", label = "Run ontology", class = "btn-success"),
                           br(),
                           fluidRow(width=12,
                                    column(width=4,
                                           br(),
                                           br(),
                                           fileInput("geneset_collection", "visit gsea-msigdb.org for your collection", accept = c("text"))),
                                    column(width = 4,
                                           br(),
                                           br(),
                                           br(),
                                           switchInput("go_plot_switch", label="switch between visualizations", onLabel="Dotplot", offLabel="Network", onStatus="primary", offStatus = "info")),
                                    column(width=4,
                                           br(),
                                           br(),
                                           br(),
                                           materialSwitch("design_plot_go", label = "change Plot design"))
                           ),
                           br(),
                           fluidRow(width=12,
                                    column(width=4,
                                           textOutput("p_correct_2")),
                                    column(width = 4,
                                           textOutput("p_cutoff_2")),
                                    column(width=4,
                                           textOutput("fc_cutoff_2"))
                           ),
                           br(),
                           conditionalPanel("!input.go_plot_switch", id="network_pane",
                                            conditionalPanel("input.design_plot_go",
                                                             visNetworkOutput("go_network_1", width = "1300px", height = "800px")
                                                             ),
                                            conditionalPanel("!input.design_plot_go",
                                                             plotOutput("go_network_2", width = "1300px", height = "800px")
                                                             )
                                            ),
                           conditionalPanel("input.go_plot_switch", id="dotplot_pane",
                                            plotlyOutput("go_dotplot", width = "1300px", height = "800px")
                                            ), 
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
  options(shiny.maxRequestSize=50*1024^2)
  

  
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  Combat_batchC <- function(i_exp_matrix, i_batch, i_model) {
    i_exp_matrix <- tryCatch(
      {
        i_exp_matrix <- ComBat(dat = i_exp_matrix,
                               batch = i_batch,
                               mod = i_model)
        return(i_exp_matrix)
      },
      error = function(cond) {
        i_exp_matrix <- ComBat(dat = i_exp_matrix,
                               batch = i_batch)
        errMsg("confounder detected! Batch and eventually your desired effect corrected. reconsider your experimental design")
        return(i_exp_matrix)
      }
    )
  }
  
  
  
  #basic error handling
  errMsg <- reactiveVal()
  output$debug <- renderPrint({
    req(errMsg())
    if (!is.null(errMsg())) {
      errMsg()
    }
  })
  
  # Help Button and content of the help menu
  observeEvent(input$help,{
    showModal(modalDialog(easyClose = T, 
                          title = "Proteomics Workbench",
                          HTML("Created by Lukas Gamp.<br> <br>
      Workflow as described in the publication: <br>
      Multiplexed single-cell proteomics using col-Median,row_Mean <br>
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
  meta_data <- eventReactive(input$update_button, {
    req(input$sample_annotation_file)
    meta_data_0 <- read.delim(input$sample_annotation_file$datapath)

    if (input$file_level == FALSE) {
      meta_data_0 <- meta_data_0[!(meta_data_0$SampleType %in%  input$selectedSampleType_to_exclude), ]
    } else {
      meta_data_0 <- meta_data_0[!(meta_data_0$Group %in%  input$selectedSampleType_to_exclude), ]
    }
    return(meta_data_0)
  })
  
  
  # read in expression set seperately for fetching raw data
  exp_set_raw <- eventReactive(input$update_button, {
    if (input$file_level == TRUE) {
      req(input$evidence_file)
      req(meta_data())
      
      meta_data_0 <- meta_data()
      evidence_data <- read.delim(input$evidence_file$datapath)
      
      #handle ID column
      rownames(evidence_data) <- evidence_data$ID
      evidence_data <- evidence_data[, !(names(evidence_data) %in% c("ID"))]
      
      #handle columns to exclude
      evidence_data <- evidence_data[ , which(meta_data_0$ID == colnames(evidence_data))]
      evidence_data[evidence_data == 0] <- NA
      return(evidence_data)
    }
  })
  
  # result of analysis pipleline
  scp <- eventReactive(input$update_button, {
    # tryCatch({
    #   errMsg(NULL)
      withProgress(message= "running analysis", value=0, {
        
        req(meta_data)
        meta_data_0 <- meta_data()
        
        if (input$file_level == FALSE) {
          req(input$evidence_file)
          evidence_data <- read.delim(input$evidence_file$datapath)
          
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
          
          
        } else {
          req(exp_set_raw())
          evidence_data <- exp_set_raw()
        } 
        
        
        incProgress(14/17, detail=paste("transforming protein data"))
        req(input$transform_base)
        transform_base_bc <- "None"
        
        if (input$transform_base == "log2") {
          print("log2 transformation")
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
          print("log10 transformation")
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
          print("sqrt transformation")
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
          print("quadratic transformation")
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
          print("boxcox transformation")
          if (input$file_level == FALSE) {
            protein_matrix <- assay(scp_0[["proteins"]])
            b <- boxcox_1(stats::lm(protein_matrix ~ 1))
            # Exact lambda
            lambda <- b$x[which.max(b$y)]
            print(paste("Lambda:", lambda))
            
            
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
            b <- boxcox_1(stats::lm(unlist(protein_matrix) ~ 1))
            # Exact lambda
            lambda <- b$x[which.max(b$y)]
            print("lambda")
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
          print("no transformation applied")
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
        if (input$norm_method == "col-Median,row_Mean" && (input$transform_base == "log2" | input$transform_base == "log10" | transform_base_bc == "log10")) {
          print("col-Median,row_Mean normalization on log transformed values")
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
            protein_matrix <- as.matrix(evidence_data)
            protein_matrix <- sweep(protein_matrix, 2, colMedians(protein_matrix), FUN="-")
            # #normalize rowwise
            protein_matrix <- sweep(protein_matrix, 1, rowMeans(protein_matrix), FUN="-")
            evidence_data <- data.frame(protein_matrix)
          }
          
        } else if (input$norm_method == "col-Median,row_Mean" && (input$transform_base != "log2" | input$transform_base != "log10" | transform_base_bc != "log10")) {
          print("col-Median,row_Mean normalization")
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
            protein_matrix <- as.matrix(evidence_data)
            protein_matrix <- sweep(protein_matrix, 2, colMedians(protein_matrix), FUN="/")
            # #normalize rowwise
            protein_matrix <- sweep(protein_matrix, 1, rowMeans(protein_matrix), FUN="/")
            evidence_data <- data.frame(protein_matrix)
          }
          
        } else if (input$norm_method == "CONSTANd") {
          print("CONSTANd normalization")
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
            evidence_data <- evidence_data$normalized_data
          }
        } else if (input$norm_method == "None") {
          print("no normalization applied")
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
        } else if (input$norm_method == "Quantile") {
          print("quantile normalization applied")
          if (input$file_level == FALSE) {
            protein_matrix <- assay(scp_0[["proteins_transf"]])
            
            sce <- getWithColData(scp_0, "proteins_transf")
            
            scp_0 <- addAssay(scp_0,
                              y = sce,
                              name = "proteins_norm")
            
            scp_0 <- addAssayLinkOneToOne(scp_0,
                                          from = "proteins_transf",
                                          to = "proteins_norm")
            
            protein_matrix <- normalizeQuantiles(protein_matrix)
            
            assay(scp_0[["proteins_norm"]]) <- protein_matrix
            
          } else {
            protein_matrix <- evidence_data
            protein_matrix <- normalizeQuantiles(evidence_data)

            evidence_data <- protein_matrix
          }
        }
        
        
        if (input$file_level == TRUE) {
          peptide_file <- unique(meta_data_0$Batch)
        }
        
        
        if (length(peptide_file) > 1) {
          incProgress(16/17, detail=paste("running missing value imputation"))
          if (input$missing_v == "KNN") {
            print("KNN missing value imputation")
            if (input$file_level == FALSE) {
              protein_matrix <- assay(scp_0[["proteins_norm"]])
              
              sce <- getWithColData(scp_0, "proteins_norm")
              
              scp_0 <- addAssay(scp_0,
                              y = sce,
                              name = "proteins_imptd")
              
              scp_0 <- addAssayLinkOneToOne(scp_0,
                                          from = "proteins_norm",
                                          to = "proteins_imptd")
              
              protein_matrix <- impute.knn(protein_matrix, 
                                           k=5, 
                                           rowmax = 1, 
                                           colmax = 1, 
                                           maxp = Inf, 
                                           rng.seed = as.numeric(gsub('[^0-9]', '', Sys.Date())))
              
              protein_matrix <- data.frame(protein_matrix$data)
              
              colnames(protein_matrix) <- colnames(assay(scp_0[["proteins_imptd"]]))
              rownames(protein_matrix) <- rownames(assay(scp_0[["proteins_imptd"]]))
              
              assay(scp_0[["proteins_imptd"]]) <- protein_matrix
            } else {
              protein_matrix <- as.matrix(evidence_data)
              protein_matrix <- impute.knn(protein_matrix,
                                           k=5,
                                           rowmax = 1,
                                           colmax = 1,
                                           maxp = Inf,
                                           rng.seed = as.numeric(gsub('[^0-9]', '', Sys.Date())))
              evidence_data <- data.frame(protein_matrix$data)
            }
          } else if (input$missing_v == "drop rows") {
            print("dropping rows with missing values")
            if (input$file_level == FALSE) {
              sce <- getWithColData(scp_0, "proteins_norm")
              
              scp_0 <- addAssay(scp_0,
                                y = sce,
                                name = "proteins_imptd")
              
              scp_0 <- addAssayLinkOneToOne(scp_0,
                                            from = "proteins_norm",
                                            to = "proteins_imptd")
              
              scp_0 <- filterNA(scp_0, pNA = 0, "proteins_imptd")
            } else {
              evidence_data <- na.omit(evidence_data)
            }
            
          } else if (input$missing_v == "replace with mean") {
            print("replacing missing values with mean")
            if (input$file_level == FALSE) {
              sce <- getWithColData(scp_0, "proteins_norm")
              
              scp_0 <- addAssay(scp_0,
                                y = sce,
                                name = "proteins_imptd")
              
              scp_0 <- addAssayLinkOneToOne(scp_0,
                                            from = "proteins_norm",
                                            to = "proteins_imptd")
              
              assay(scp_0[["proteins_imptd"]]) <- replace(assay(scp_0[["proteins_imptd"]]), is.na(assay(scp_0[["proteins_imptd"]])), mean(assay(scp_0[["proteins_imptd"]]), na.rm = TRUE))
            } else {
              evidence_data <- replace(evidence_data, is.na(evidence_data), mean(unlist(evidence_data), na.rm = TRUE))
            }
            
          } else if (input$missing_v == "replace with median") {
            print("replacing missing values with median")
            if (input$file_level == FALSE) {
              sce <- getWithColData(scp_0, "proteins_norm")
              
              scp_0 <- addAssay(scp_0,
                                y = sce,
                                name = "proteins_imptd")
              
              scp_0 <- addAssayLinkOneToOne(scp_0,
                                            from = "proteins_norm",
                                            to = "proteins_imptd")
              
              assay(scp_0[["proteins_imptd"]]) <- replace(assay(scp_0[["proteins_imptd"]]), is.na(assay(scp_0[["proteins_imptd"]])), median(assay(scp_0[["proteins_imptd"]]), na.rm = TRUE))
            } else {
              evidence_data <- replace(evidence_data, is.na(evidence_data), median(unlist(evidence_data), na.rm = TRUE))
            }
          }
          
          incProgress(16/17, detail=paste("running batch correction"))
          if (input$batch_c == "ComBat") {
            print("running ComBat")
            if (input$file_level == FALSE) {
              sce <- getWithColData(scp_0, "proteins_imptd")
              batch <- colData(sce)$Raw.file
              model <- model.matrix(~SampleType, data = colData(sce))
              
              assay(sce) <- Combat_batchC(assay(sce), batch, model) 
              
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
            } else {
              batch <- meta_data_0$Batch
              model <- model.matrix(~Group, data=meta_data_0)
              
              evidence_data <- Combat_batchC(evidence_data, batch, model)
            }
          } else if (input$batch_c == "None") {
            print("no batch correction applied")
            if (input$file_level == FALSE) {
              sce <- getWithColData(scp_0, "proteins_imptd")
              
              scp_0 <- addAssay(scp_0,
                                y = sce,
                                name = "proteins_dim_red")
              
              scp_0 <- addAssayLinkOneToOne(scp_0,
                                            from = "proteins_imptd",
                                            to = "proteins_dim_red")          
            } else {
              evidence_data <- evidence_data
            }
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
        if (input$file_level == FALSE) {
          scp_0[["proteins_dim_red"]] <- scater::runPCA(scp_0[["proteins_dim_red"]],
                                                        ncomponents = 3,
                                                        ntop = Inf,
                                                        scale = TRUE,
                                                        exprs_values = 1,
                                                        name = "PCA")
          
          scp_0[["proteins_dim_red"]] <- runUMAP(scp_0[["proteins_dim_red"]],
                                                 ncomponents = 3,
                                                 ntop = Inf,
                                                 scale = TRUE,
                                                 exprs_values = 1,
                                                 n_neighbors = 3,
                                                 dimred = "PCA",
                                                 name = "UMAP")
        }
        incProgress(17/17, detail=paste("analysis finish"))
      })
      if (input$file_level == FALSE) {
        return(scp_0)
      } else {
        return(evidence_data)
      }
    # }, error = function(err) {
    #   print("error handler")
    #   errMsg("Whoops something went wrong")
    # },
    # finally = invalidateLater(1))
  })
  
  
  # expression matrix is used by plot functions
  exp_matrix <- reactive({
    req(scp())
    if (input$file_level == FALSE) {
      scp_0 <- scp()
      exp_matrix_0 <- data.frame(assay(scp_0[["proteins_dim_red"]]))
      colnames(exp_matrix_0) <- scp_0$SampleType  
    } else {
      req(meta_data())
      exp_matrix_0 <- data.frame(scp())
      meta_data_0 <- meta_data()
      colnames(exp_matrix_0) <- meta_data_0$Group
    }
    return(exp_matrix_0)
  })
  
  # trigger for the CONSTANd dependency MA plot
  CONSTANd_trigger <- reactive(list(input$update_button, input$norm_method))
  ### graphical outputs and user interface
  
  
  # observer for protein list
  observe({
    req(protein_list())
    updateSelectizeInput(session, "selectedProtein", choices = protein_list())
  })
  
  # reactive element for the protein list --> update if the scp object changes
  protein_list <- eventReactive(input$update_button, {
    req(scp())
    if (input$file_level == FALSE) {
      req(scp()[["proteins"]])
      scp_0 <- scp()
      list <- rowData(scp_0)[["proteins"]][,1]
    } else {
      list <- rownames(scp())
    }
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
  comp_list <- eventReactive(input$update_button, {
    if (input$file_level == FALSE) {
      req(scp())
      scp_0 <- scp()
      list <- apply(combn(unique(scp_0$SampleType),2),2,paste, collapse="-")
    } else {
      req(meta_data())
      list <- apply(combn(unique(meta_data()$Group),2),2,paste, collapse="-")
    }
    return(list)
  })
  
  columns <- eventReactive(input$update_button, {
    if (input$file_level == FALSE) {
      req(scp())
      scp_0 <- scp()
      list <- colnames(colData(scp_0))
    } else {
      req(meta_data())
      list <- colnames(meta_data())
    }
    return(list)
  })
  
  # observer for sample type exclusion
  observe({
    req(sample_types())
    updateSelectInput(session, "selectedSampleType_to_exclude", choices = sample_types())
  })
  
  sample_types <- eventReactive(input$update_button, {
    if (input$file_level == FALSE) {
      req(scp())
      scp_0 <- scp()
      list <- scp_0$SampleType
    } else {
      req(meta_data())
      meta_data_0 <- meta_data()
      list <- meta_data_0$Group
    }
    return(list)
  })
  
  ## Output of analysis pipeline 
  # pathway of the data first tab
  # include pipeline for expression set
  output$overview_plot <- renderPlot({
    if (!is.null(scp()) && input$file_level == FALSE) {
      plot(scp()) }
  })
  
  # summary barchart second tab
  output$summary_bar <- renderPlot({
    if (input$file_level == FALSE) {
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
    } else {
      barplot(log10(dim(scp())[1]), main = "log Count of rows", horiz=T, )
    }
  })
  
  # Reporter Ion intensity visualisation
  #observer for color_variable
  observe({
    req(columns())
    updateSelectInput(session, "color_variable_ri", choices=columns())
  })  
  # Boxplot of reporter ion intensity third tab
  output$RI_intensity <- renderPlot({
    
    if (input$file_level == FALSE) {
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
    } else{
      meta_data_0 <- meta_data()
      evidence_data <- scp()
      
      evidence_data_medians <- data.frame(median_intensity = colMedians(as.matrix(evidence_data)))
      evidence_data_medians %>%
        ggplot() +
        aes(x = median_intensity, 
            y = meta_data_0[[as.character(input$color_variable_ri)]],
            fill = meta_data_0[[as.character(input$color_variable_ri)]]) +
        geom_boxplot() +
        scale_x_log10() +
        labs(fill=as.character(input$color_variable_ri), y=as.character(input$color_variable_ri))
    }
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
    
    if (input$file_level == FALSE) {
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
          labs(fill=as.character(input$color_variable_cv), y=as.character(input$color_variable_cv), title = "Covariance over razor proteins")  
      } else {
        getWithColData(scp_0, peptide_file) %>%
          colData %>%
          data.frame %>%
          ggplot(aes(x = MedianCV,
                     fill = get(input$color_variable_cv))) +
          geom_boxplot()+
          labs(fill=as.character(input$color_variable_cv), y=as.character(input$color_variable_cv), title = "Covariance over razor proteins")
      }
    }
  })
  
  output$corr_matrix <- renderPlot({
    scp_0 <- scp()
    if (input$file_level == FALSE) {
      heatmap(cor(t(assay(scp_0[["proteins_dim_red"]]))))
    } else {
      heatmap(cor(t(scp_0)))
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
  
  
  dimred_pca_expr_set <- eventReactive(input$update_button, {
    if (input$file_level == TRUE) {
      req(scp())
      scp_0 <- scp()
      dimred_pca <- calculatePCA(scp_0,                                          
                                 ncomponents = 5,
                                 ntop = Inf,
                                 scale = TRUE)
      
      dimred_pca <- data.frame(dimred_pca)
      return(dimred_pca)
    }
  })
  
  
  # principle component analysis in fifth tab
  output$PCA <- renderPlot({
    
    extract_null <- function(variable) if (variable == "NULL") {
      return(NULL)
    } else {
      if (input$file_level == FALSE) {
        return(variable)
      } else {
        return(meta_data_0[[as.character(variable)]])
      }
    }
    meta_data_0 <- meta_data()
    scp_0 <- scp()
    
    color_variable_dim_red <- extract_null(input$color_variable_dim_red)
    shape_variable_dim_red <- extract_null(input$shape_variable_dim_red)
    size_variable_dim_red <- extract_null(input$size_variable_dim_red)
    
    if (input$file_level == FALSE) {
      plotReducedDim(scp_0[["proteins_dim_red"]],
                     dimred = "PCA",
                     colour_by = color_variable_dim_red,
                     shape_by = shape_variable_dim_red,
                     size_by = size_variable_dim_red,
                     point_alpha = 1,
                     point_size=3)
    } else {
      
      dimred_pca <- dimred_pca_expr_set()
      
      ggplot(dimred_pca, aes(x=PC1, 
                             y=PC2, 
                             color= color_variable_dim_red,
                             shape = shape_variable_dim_red,
                             size = size_variable_dim_red)) +
        geom_point(alpha=3/4)
    }
  })
  
  output$PCA_third_dim <- renderPlotly({
    
    extract_null <- function(variable) if (variable == "NULL") {
      return(NULL)
    } else {
      if (input$file_level == FALSE) {
        return(scp_0[["proteins_dim_red"]][[variable]])
      } else {
        return(meta_data_0[[as.character(variable)]])
      }
    }
    
    meta_data_0 <- meta_data()
    scp_0 <- scp()
    
    
    color_variable_dim_red <- extract_null(input$color_variable_dim_red)
    shape_variable_dim_red <- extract_null(input$shape_variable_dim_red)
    size_variable_dim_red <- extract_null(input$size_variable_dim_red)
    
    
    if (input$file_level == FALSE) {
      
      dimred_pca <- data.frame(reducedDim(scp_0[["proteins_dim_red"]], "PCA"))
      plot_ly(x=dimred_pca$PC1,
              y=dimred_pca$PC2,
              z=dimred_pca$PC3,
              type="scatter3d",
              mode="markers",
              color=color_variable_dim_red)
      
    } else {
      
      dimred_pca <- dimred_pca_expr_set()
      
      plot_ly(x=dimred_pca$PC1,
              y=dimred_pca$PC2,
              z=dimred_pca$PC3,
              type="scatter3d",
              mode="markers",
              color=color_variable_dim_red)
      
    }
  })
  
  dimred_umap_expr_set <- eventReactive(input$update_button, {
    if (input$file_level == TRUE) {
      req(scp())
      scp_0 <- scp()
      dimred_umap <- calculateUMAP(scp_0,
                                   ncomponents = 3,
                                   ntop = Inf,
                                   scale = TRUE)
      
      dimred_umap <- data.frame(dimred_umap)
      return(dimred_umap)
    }
  })
  
  
  
  # Umap dimensionality reduction in fith tab
  output$UMAP <- renderPlot({
    
    extract_null <- function(variable) if (variable == "NULL") {
      return(NULL)
    } else {
      if (input$file_level == FALSE) {
        return(variable)
      } else {
        return(meta_data_0[[as.character(variable)]])
      }
    }
    meta_data_0 <- meta_data()
    scp_0 <- scp()
    
    color_variable_dim_red <- extract_null(input$color_variable_dim_red)
    shape_variable_dim_red <- extract_null(input$shape_variable_dim_red)
    size_variable_dim_red <- extract_null(input$size_variable_dim_red)
    
    if (input$file_level == FALSE) {
      plotReducedDim(scp_0[["proteins_dim_red"]],
                     dimred = "UMAP",
                     colour_by = color_variable_dim_red,
                     shape_by = shape_variable_dim_red,
                     size_by = size_variable_dim_red,
                     point_alpha = 1,
                     point_size = 3)      
    } else {
      dimred_umap <- dimred_umap_expr_set()
      
      ggplot(dimred_umap, aes(x=UMAP1, 
                              y=UMAP2, 
                              color= color_variable_dim_red,
                              shape = shape_variable_dim_red,
                              size = size_variable_dim_red)) +
        geom_point(alpha=3/4)
    } 
  })
  
  output$UMAP_third_dim <-renderPlotly({
    
    extract_null <- function(variable) if (variable == "NULL") {
      return(NULL)
    } else {
      if (input$file_level == FALSE) {
        return(scp_0[["proteins_dim_red"]][[variable]])
      } else {
        return(meta_data_0[[as.character(variable)]])
      }
    }
    
    meta_data_0 <- meta_data()
    scp_0 <- scp()    
    
    color_variable_dim_red <- extract_null(input$color_variable_dim_red)
    shape_variable_dim_red <- extract_null(input$shape_variable_dim_red)
    size_variable_dim_red <- extract_null(input$size_variable_dim_red)
    
    
    
    if (input$file_level == FALSE) {
      dimred_umap <- data.frame(reducedDim(scp_0[["proteins_dim_red"]], "UMAP"))
      plot_ly(x=dimred_umap$UMAP1,
              y=dimred_umap$UMAP2,
              z=dimred_umap$UMAP3,
              type="scatter3d",
              mode="markers",
              color=color_variable_dim_red)  
    } else {
      dimred_umap<- dimred_umap_expr_set() 
      
      plot_ly(x=dimred_umap$UMAP1,
              y=dimred_umap$UMAP2,
              z=dimred_umap$UMAP3,
              type="scatter3d",
              mode="markers",
              color=color_variable_dim_red)
      
    } 
  })
  
  observe({
    req(columns())
    updateSelectInput(session, "color_variable_feature", choices=columns())
  })  
  
  # protein wise visualisation 
  output$feature_subset <- renderPlot({
    scp_0 <- scp()
    
    
    channelstring <- gsub("[0-9]{1,2}$","", scp_0$Channel[1]) 
    
    if (input$file_level == FALSE) {
      subsetByFeature(scp_0, input$selectedProtein) %>%
        ## Format the `QFeatures` to a long format table
        longFormat(colvars = c("Raw.file", "SampleType", "Channel")) %>%
        data.frame %>%
        ## This is used to preserve ordering of the samples and assays in ggplot2
        mutate(assay = factor(assay, levels = names(scp_0)),
               Channel = sub(channelstring, "", Channel)) %>%
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
    } else {
      extract_null <- function(variable) if (variable == "NULL") {
        return(NULL)
      } else {
        if (input$file_level == FALSE) {
          return(scp_0[["proteins_dim_red"]][[variable]])
        } else {
          return(meta_data_0[[as.character(variable)]])
        }
      }
      
      
      req(meta_data())
      meta_data_0 <- meta_data()
      
      to_plot <- data.frame(t(scp_0[input$selectedProtein, ]))
      
      color_variable_feature <- extract_null(input$color_variable_feature)
      
      d = data.frame(melt(to_plot))
      
      ggplot(data = d,
             mapping = aes(x = variable, y = value, fill=color_variable_feature)) + 
        geom_col(position = position_dodge()) +
        labs(y=paste("intensity of", input$selectedProtein))
    }
  })
  
  # plot MA plot for the constand method
  output$MAplot <- renderPlot({
    req(input$selectedComp)
    user_choice <- input$selectedComp
    
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
      
      use <- is.finite(A)
      A <- A[use]
      M <- M[use]
      
      # plot
      print(var(M))
      plot(A, M, xlab=xlab, cex.lab=cex, cex.axis=cex)
      # rolling average
      if (showavg) { lines(lowess(M~A), col='red', lwd=5) }
    }
    
    if (input$file_level == FALSE) {
      req(scp())
      scp_0 <- scp() 
      
      scp_0 <- filterNA(scp_0, pNA = 0, "proteins_transf")
      
      
      # split user choice of comp back to sample types
      user_choice_vector <- strsplit(user_choice, split = "-")
      # and assign them to a variable
      choice_A <- user_choice_vector[[1]][1]
      choice_B <- user_choice_vector[[1]][2]
      
      #find the row indeces of corresponding to the individual sample types
      st_indeces <- split(seq_along(scp_0$SampleType), scp_0$SampleType)
      index_A <- st_indeces[choice_A]
      index_B <- st_indeces[choice_B]
      
      MAplot(assay(scp_0[["proteins_transf"]][,index_A[[1]]]), assay(scp_0[["proteins_transf"]][,index_B[[1]]]))
    } else {
      req(exp_matrix())
      req(meta_data())
      exp_matrix_0 <- as.matrix(exp_matrix())
      meta_data_0 <- meta_data()
      
      # split user choice of comp back to sample types
      user_choice_vector <- strsplit(user_choice, split = "-")
      # and assign them to a variable
      choice_A <- user_choice_vector[[1]][1]
      choice_B <- user_choice_vector[[1]][2]
      
      #find the row indeces of corresponding to the individual sample types
      st_indeces <- split(seq_along(meta_data_0$Group), meta_data_0$Group)
      index_A <- st_indeces[choice_A]
      index_B <- st_indeces[choice_B]
      
      MAplot(exp_matrix_0[,index_A[[1]]], exp_matrix_0[,index_B[[1]]])
    }
    
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
  
  
  observe({
    req(sample_types())
    updateSelectInput(session, "selectedSampleType_normval", choices = sample_types())
  })
  
  norm_method_val <- eventReactive(input$update_button,{
    req(input$norm_method)
    return(paste(input$norm_method))
  })
  
  transform_base_val <- eventReactive(input$update_button,{
    req(input$transform_base)
    return(paste(input$transform_base))
  })
  
  
  ## Validation of the procedures
  output$norm_val_plot_intra_group <- renderPlot({
    
    req(scp())
    scp_0 <- scp()
    
    if (input$file_level == FALSE) {
      # select data according to the procedure undertaken
      data_final <- assay(scp_0[["proteins_dim_red"]])
      data_raw <- assay(scp_0[["proteins"]])
      
      # select data according to groups to observe mutual information within
      data_final <- data_final[, which(scp_0$SampleType == input$selectedSampleType_normval)]
      data_raw <- data_raw[, which(scp_0$SampleType == input$selectedSampleType_normval)]
      
      req(input$selectedSampleType_normval)
      req(norm_method_val())
      req(transform_base_val())
      
    } else {
      data_final <- scp_0
      
      req(exp_set_raw())
      data_raw <- exp_set_raw()
      
      req(meta_data())
      meta_data_0 <- meta_data()
      
      # select data according to groups to observe mutual information within
      data_final <- data_final[, which(meta_data_0$Group == input$selectedSampleType_normval)]
      data_raw <- data_raw[, which(meta_data_0$Group == input$selectedSampleType_normval)]
      
    }
    
    
    
    # discretize by equal frequencies
    data_final <- discretize(data_final)
    data_raw <- discretize(data_raw)
    
    # mutual information is returned in nats
    mi_final <- mutinformation(data_final)
    mi_raw <- mutinformation(data_raw)
    
    mi_final_stack <- data.frame(stack(mi_final))
    mi_raw_stack <- data.frame(stack(mi_raw))
    
    mi_total <- cbind(mi_final_stack$value, mi_raw_stack$value)
    colnames(mi_total) <- c("MI_final", "MI_raw")
    
    to_plot <- melt(mi_total)
    
    ggplot(to_plot, aes(x=Var2, y=value)) +
      geom_boxplot(aes(fill=Var2)) +
      ggtitle(paste("Change in mutual information after transformation:", 
                    transform_base_val(), 
                    "and norm_method:", norm_method_val())) +
      xlab(paste(input$selectedSampleType_normval)) +
      ylab("natural unit of information (nat)")+ 
      guides(fill=guide_legend(title=""))
  })
  
  
  observe({
    updateSelectInput(session, "selectedComp_normval", choices = comp_list())
  })
  
  
  output$norm_val_plot_inter_group <- renderPlot({
    
    req(scp())
    scp_0 <- scp()
    
    req(input$selectedComp_normval)
    user_choice <- input$selectedComp_normval
    
    # split user choice of comp back to sample types
    user_choice_vector <- strsplit(user_choice, split = "-")
    # and assign them to a variable
    choice_A <- user_choice_vector[[1]][1]
    choice_B <- user_choice_vector[[1]][2]
    
    if (input$file_level == FALSE) {
      
      # select data according to the procedure undertaken
      data_final <- assay(scp_0[["proteins_dim_red"]])
      data_raw <- assay(scp_0[["proteins"]])
      
      # select data according to groups to observe mutual information within
      data_final_group1 <- discretize(data_final[, which(scp_0$SampleType == choice_A)])
      data_raw_group1 <- discretize(data_raw[, which(scp_0$SampleType == choice_A)])
      
      data_final_group2 <- discretize(data_final[, which(scp_0$SampleType == choice_B)])
      data_raw_group2 <- discretize(data_raw[, which(scp_0$SampleType == choice_B)])
      
    } else {
      data_final <- scp_0
      
      req(exp_set_raw())
      data_raw <- exp_set_raw()
      
      req(meta_data())
      meta_data_0 <- meta_data()
      
      # select data according to groups to observe mutual information within
      data_final_group1 <- discretize(data_final[, which(meta_data_0$Group == choice_A)])
      data_raw_group1 <- discretize(data_raw[, which(meta_data_0$Group == choice_A)])
      
      data_final_group2 <- discretize(data_final[, which(meta_data_0$Group == choice_B)])
      data_raw_group2 <- discretize(data_raw[, which(meta_data_0$Group == choice_B)])
      
    }
    
    data_final <- cbind(data_final_group1, data_final_group2)
    data_raw <- cbind(data_raw_group1, data_raw_group2)
    
    mi_final <- mutinformation(data_final)
    mi_raw <- mutinformation(data_raw)
    
    extr_val_above_diag <- function(mat) {
      vec <- c()
      for(row in 1:nrow(mat))
      {
        
        # looping over columns
        for(col in 1:ncol(mat))
        {
          # if column number is greater than row
          if(col > row)
          {
            # printing the element of the matrix
            vec <- append(vec, mat[row,col])
          }
        }
      }
      return(vec)
    }
    
    mi_final <- extr_val_above_diag(mi_final)
    mi_raw <- extr_val_above_diag(mi_raw)
    
    mi_total <- data.frame(MI_final=mi_final, MI_raw=mi_raw)
    
    to_plot <- melt(mi_total)
    
    req(norm_method_val())
    req(transform_base_val())
    
    ggplot(to_plot, aes(x=variable, y=value)) +
      geom_boxplot(aes(fill=variable)) +
      ggtitle(paste("Change in mutual information after transformation:", 
                    transform_base_val(), 
                    "and norm_method:", norm_method_val())) +
      xlab(paste(input$selectedComp_normval)) +
      ylab("natural unit of information (nat)")+ 
      guides(fill=guide_legend(title=""))
  })
  
  
  
  
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
      meta_data_0 <- meta_data()
      
      
      incProgress(2/3, detail=paste("creating linear model"))
      req(input$model_design)
      # Create a design matrix
      if (input$model_design == "All pairwise comparison") {
        if (input$file_level == FALSE) {
          SampleType <- scp_0$SampleType
          design <- model.matrix(~0+SampleType)
          colnames(design) <- sub("SampleType", "", colnames(design))
          fit <- lmFit(exp_matrix_0, design)
        } else {
          Group <- meta_data_0$Group
          design <- model.matrix(~0+Group)
          colnames(design) <- sub("Group", "", colnames(design))
          fit <- lmFit(exp_matrix_0, design)
        }
        
      }
      #Differential Expression with defined Contrasts
      else if (input$model_design == "Differential Expression with defined Contrasts") {
        if (input$file_level == FALSE) {
          # fetch user selection
          req(input$selectedComp_stat)
          
          SampleType <- scp_0$SampleType
          design <- model.matrix(~0+SampleType)
          colnames(design) <- sub("SampleType", "", colnames(design))
          
          fit <- lmFit(exp_matrix_0, design)
          
          user_contrast <- paste(input$selectedComp_stat, collapse = "-")
          cont_matrix <- makeContrasts(contrasts=user_contrast,levels=colnames(design))
          
          fit <- contrasts.fit(fit, cont_matrix)
        } else {
          req(input$selectedComp_stat)
          
          Group <- meta_data_0$Group
          design <- model.matrix(~0+Group)
          colnames(design) <- sub("Group", "", colnames(design))
          
          fit <- lmFit(exp_matrix_0, design)
          
          user_contrast <- paste(input$selectedComp_stat, collapse = "-")
          cont_matrix <- makeContrasts(contrasts=user_contrast,levels=colnames(design))
          
          fit <- contrasts.fit(fit, cont_matrix)
        }
      }
      else if (input$model_design == "Multi factor additivity") {
        if (input$file_level == FALSE) {
          req(input$col_factors)
          req(input$selectedComp_stat)
          
          fetched_factor <- colData(scp_0)[input$col_factors]
          design_frame <- cbind(scp_0$SampleType, fetched_factor)
          
          colnames(design_frame) <- c("SampleType", input$col_factors)
          
          
          design <- model.matrix(~0+ . , data=design_frame)
          colnames(design)[1:length(unique(design_frame$SampleType))] <- sub("SampleType", "", colnames(design)[1:length(unique(design_frame$SampleType))])
          
          fit <- lmFit(exp_matrix_0, design)
          
          user_contrast <- paste(input$selectedComp_stat, collapse = "-")
          cont_matrix <- makeContrasts(contrasts=user_contrast,levels=colnames(design))
          
          fit <- contrasts.fit(fit, cont_matrix)
          
        } else {
          req(input$col_factors)
          req(input$selectedComp_stat)
          
          fetched_factor <- meta_data_0[input$col_factors]
          design_frame <- cbind(meta_data_0$Group, fetched_factor)
          
          colnames(design_frame) <- c("Group", input$col_factors)
          
          design <- model.matrix(~0+ . , data=design_frame)
          colnames(design)[1:length(unique(design_frame$Group))] <- sub("Group", "", colnames(design)[1:length(unique(design_frame$Group))])
          
          fit <- lmFit(exp_matrix_0, design)
 
          user_contrast <- paste(input$selectedComp_stat, collapse = "-")
          cont_matrix <- makeContrasts(contrasts=user_contrast,levels=colnames(design))

          fit <- contrasts.fit(fit, cont_matrix)
          
        }
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
    selectInput("selectedComp_stat", "choose your contrast of interest", "", multiple = T)
  })
  
  # Ui for multifactorial model
  output$add_factor <- renderUI({
    req(input$model_design == "Multi factor additivity")
    selectInput("col_factors", "choose additional factor(s)", "", multiple = T)
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
    
    toptable$Significance <- "not significant"
    toptable$Significance[(toptable$logFC > FCcutoff)] <- "upregulated"
    toptable$Significance[(toptable$logFC < -1*(FCcutoff))] <- "downregulated"
    toptable$Significance[(toptable$adj.P.Val < pCutoff)] <- "p-value"
    
    toptable$Significance[(toptable$adj.P.Val < pCutoff) & (toptable$logFC > FCcutoff)] <- "significantly upregulated"
    
    toptable$Significance[(toptable$adj.P.Val < pCutoff) & (toptable$logFC < -1*(FCcutoff))] <- "significantly downregulated"
    
    toptable$Significance <- factor(toptable$Sig, levels = c("not significant",
                                                             "upregulated",
                                                             "downregulated",
                                                             "p-value", 
                                                             "significantly upregulated",
                                                             "significantly downregulated"
    ))
    
    req(input$p_value_correction)
    if (input$p_value_correction == "none") {
      p_val_correct <- "p-value"
    } else {
      p_val_correct <- "adj.p.value"
    }
    
    p <- ggplot(toptable, aes(x = logFC, y = -log10(adj.P.Val), text=protein)) +
      geom_point(aes(color = Significance), alpha = 3/4, shape = 19, size = 1.5, na.rm = TRUE) +
      theme_bw() +
      theme(legend.position = "top") +
      scale_color_manual(values = c("not significant" = "grey",
                                    "upregulated" = "#44AA99",
                                    "downregulated" = "#AA4499",
                                    "p-value" = "#88CCEE", 
                                    "significantly upregulated" = "green",
                                    "significantly downregulated" = "red")) +
      labs(y=p_val_correct)
    
    ggplotly(p, tooltip = "text")
    
    
  })
  
  # reactive element for toptable
  protein_table <- reactive({
    req(stat_result())
    req(input$chosen_coef)
    data.frame(topTable(stat_result(), number = Inf, adjust.method = input$p_value_correction, coef = input$chosen_coef))
  })
  
  # show significant proteins in the table
  output$protein_table <- DT::renderDataTable({
    req(protein_table())
    protein_table()# %>%
    #      mutate_if(is.numeric, round, digits=10)
  })
  
  # venn diagram for significant proteins
  output$venn_diagram <- renderPlot({
    req(stat_result())
    vennDiagram(decideTests(stat_result(), p.value=input$p_value_cutoff, adjust.method=input$p_value_correction))
  })
  
  
  # gene set enrichment analysis 
  
  result_kegg <- eventReactive(input$run_kegg, {
    withProgress(message = "running pathway analysis", value = 0, {
      incProgress(1/3, detail = paste("reading data"))
      req(protein_table())
      tt <- protein_table()
      
      # filter the data according to user selection
      req(input$p_value_cutoff)
      req(input$fold_change_cutoff)
      p_cut <- input$p_value_cutoff
      fc_cut <- input$fold_change_cutoff
      
      mask <- tt$adj.P.Val < p_cut & 
        abs(tt$logFC) > fc_cut
      
      deGenes <- rownames(tt)[mask]
      
      if (input$panel_as_universe) {
        background <- rownames(tt)
      } else {
        background <- NULL
      }
      
      if (!is.null(input$custom_universe) & !(input$panel_as_universe)) {
        background <- read.delim(input$custom_universe$datapath)
        background <- as.vector(background[1])
        background <- background[[1]]
      } else {
        background <- background
      }
      
      req(input$p_value_correction)
      p_correct <- input$p_value_correction
      
      incProgress(2/3, detail = paste("running pathway enrichment"))
      ans.kegg <- enrichKEGG(
        deGenes,
        organism = "hsa",
        keyType = "uniprot",
        pvalueCutoff = p_cut,
        pAdjustMethod = p_correct,
        minGSSize = 10,
        maxGSSize = Inf,
        qvalueCutoff = 0.2,
        use_internal_data = FALSE,
        universe = background
      )
      incProgress(3/3, detail = paste("success"))
    })
    return(ans.kegg)
  })
  
  output$p_correct_1 <- output$p_correct_2 <- renderText({
    req(input$p_value_correction)
    paste(
      "p-value correction method:",
      input$p_value_correction
    )
  })
  
  
  output$p_cutoff_1 <- output$p_cutoff_2  <- renderText({
    req(input$p_value_cutoff)
    paste(
      "p-value cutoff:",
      input$p_value_cutoff
    )
    
  })
  output$fc_cutoff_1 <- output$fc_cutoff_2 <- renderText({
    req(input$fold_change_cutoff)
    paste(
      "fold change cutoff:",
      input$fold_change_cutoff
    )
  })
  
  output$gsea<- renderPlotly({
    req(result_kegg())
    ans.kegg <- result_kegg()
    tab.kegg <- as.data.frame(ans.kegg)
    
    if (input$design_plot_gsea) {
      
      tab.kegg$denominator <- as.numeric(gsub("^\\d+/(\\d+)$", "\\1", tab.kegg$GeneRatio))
      tab.kegg$decimal_gene_ratio <- tab.kegg$Count / tab.kegg$denominator
      
      plot_ly(data=tab.kegg,
              x=~decimal_gene_ratio,
              y=~Description,
              type = "scatter",
              color= ~p.adjust,
              size=~decimal_gene_ratio,
              text=~GeneRatio,
              hovertemplate= paste('%{y}',
                                   '<br>Gene ratio: %{text}<br>
                                   <extra></extra>')
      ) %>%
        layout(xaxis=list(
          title="Gene Ratio",
          range=c(0,1)
        ))
    } else {
      graphics::barplot(ans.kegg, showCategory=10)
    }
  })
  
  
  gsea_network_data <- reactive({
    req(result_kegg())
    ans.kegg <- result_kegg()
    
    nw <- cnetplot(ans.kegg)
    return(nw)
  })
  
  
  output$gsea_network_1 <- renderVisNetwork({
    if (input$design_plot_gsea & !input$pw_plot_switch) {
      withProgress(message = "Plotting network", value=0,{
        incProgress(1/3, detail=paste("reading data"))
        req(gsea_network_data())
        nw_data <- gsea_network_data()
        
        graph_info <- attributes(nw_data$data)$graph
        incProgress(2/3, detail = paste("rendering network"))
        nw <- visNetwork::visIgraph(graph_info, physics = T, smooth = T)
        
        incProgress(3/3, detail = paste("success"))
      })
      nw 
    }
  })
  
  output$gsea_network_2 <- renderPlot({
    withProgress(message = "Plotting network", value=0,{
      incProgress(1/3, detail=paste("reading data"))
      req(gsea_network_data())
      nw_data <- gsea_network_data()
      
      incProgress(2/3, detail = paste("rendering network"))
      nw <- nw_data
      
      incProgress(3/3, detail = paste("success"))
    })
    nw
  })
  
  result_piano <- eventReactive(input$go_ontology, {
    withProgress(message = "running protein set analysis ontology", value = 0, {
      req(protein_table())
      tt <- protein_table()
      
      # parse data
      myPval <- data.frame(p.val=tt$adj.P.Val, t.val=tt$t, logFC=tt$logFC)
      myPval$UNIPROT <- rownames(tt)

      incProgress(1/4, detail = paste("fetching entrez IDs"))
      entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, myPval$UNIPROT, "ENTREZID", "UNIPROT")
      myPval <- merge(myPval, entrez_ids, by="UNIPROT")
      
      # drop duplicates and nas
      myPval <- myPval[which(!is.na(myPval$ENTREZID)) , ]
      myPval <- myPval[!duplicated(myPval$ENTREZID) , ]
      
      rownames(myPval) <- myPval$ENTREZID
      
      incProgress(2/4, detail= paste("loading protein set collection"))
      req(input$geneset_collection)
      myGSC <- loadGSC(input$geneset_collection$datapath)
      
      
      incProgress(3/4, detail=paste("running protein set analysis"))
      logFCs <- data.frame(myPval$logFC)
      rownames(logFCs) <- rownames(myPval)
      
      pVals <- data.frame(myPval$p.val)
      rownames(pVals) <- rownames(myPval)
      
      req(input$p_value_correction)
      p_correct <- input$p_value_correction
      
     cores <- detectCores()

     gsaRes <- runGSA(geneLevelStats = pVals,
                      directions = logFCs,
                      gsc=myGSC,
                      adjMethod = p_correct,
                      geneSetStat = "fisher",
                      ncpus = cores
                      )
      incProgress(4/4, detail=paste("success"))
      
    })
    return(gsaRes)
  })
  
  
  
  output$go_network_1 <- renderVisNetwork({
    if (input$design_plot_go & !input$go_plot_switch) {
      withProgress(message = "Plotting network", value=0,{
        incProgress(1/3, detail=paste("reading data"))
        req(result_piano())
        gsaRes <- result_piano()
        
        req(input$p_value_cutoff)
        p_cut <- input$p_value_cutoff
        
        incProgress(2/3, detail = paste("rendering network"))
        nw <- networkPlot2(gsaRes, class="non", significance = p_cut, shiny=T)
        
        incProgress(3/3, detail = paste("success"))
      })
      nw
    }
  })
  
  output$go_network_2 <- renderPlot({
      withProgress(message = "Plotting network", value=0,{
        incProgress(1/3, detail=paste("reading data"))
        req(result_piano())
        gsaRes <- result_piano()
        
        req(input$p_value_cutoff)
        p_cut <- input$p_value_cutoff
        
        incProgress(2/3, detail = paste("rendering network"))
        nw <- networkPlot(gsaRes, class="non", significance = p_cut)
        
        incProgress(3/3, detail = paste("success"))
      })
      nw 
  })
  
  
  
  output$go_dotplot <- renderPlotly({
    withProgress(message = "Plotting dotplot", value=0, 
                 {
                   incProgress(1/3, detail = paste("reading data"))
                   
                   req(result_piano())
                   gsaRes <- result_piano()
                   
                   req(input$p_value_cutoff)
                   p_cut <- input$p_value_cutoff
                   
                   
                   gsa_results <- GSAsummaryTable(gsaRes = gsaRes)
                   gsa_results$name_wo_suff <- gsub("_UP*", "", gsa_results$Name)
                   gsa_results$name_wo_suff <- gsub("_DN*", "", gsa_results$name_wo_suff)
                   
                   
                   gsa_results_total_count <- data.frame(aggregate(gsa_results$`Genes (tot)`, list(gsa_results$name_wo_suff), sum))
                   colnames(gsa_results_total_count) <- c("Description", "count")
                   
                   
                   gsa_results_total_pvals <- data.frame(aggregate(gsa_results$`p adj (non-dir.)`, list(gsa_results$name_wo_suff), mean))
                   colnames(gsa_results_total_pvals) <- c("Description", "p.adjust")
                   
                   
                   gsa_results_total <- cbind(gsa_results_total_count, gsa_results_total_pvals)
                   
                   
                   gsa_results_total$ratio <- paste(as.character(gsa_results_total$count), "/" ,as.character(length(gsa_results_total$Description)))

                                      
                   gsa_results_total$decimal_ratio <- gsa_results_total$count / length(gsa_results_total$Description)
                   
                   
                   gsa_results_total <- gsa_results_total[, !duplicated(colnames(gsa_results_total))]
                   
                   mask <- gsa_results_total$p.adjust < p_cut

                   gsa_results_total <- gsa_results_total[mask, ]



                   if (input$design_plot_go) {
                     incProgress(2/3, detail=paste("rendering dotplot"))
                     p1 <- plot_ly(data=gsa_results_total,
                             x=~count,
                             y=~Description,
                             type = "scatter",
                             color= ~p.adjust,
                             size=~decimal_ratio,
                             text=~ratio,
                             hovertemplate= paste('%{y}', '<br>Gene ratio: %{text}<br><extra></extra>')) %>%
                       layout(xaxis=list(
                         title="Gene Ratio"
                       ))
                   } else {
                     incProgress(2/3, detail=paste("rendering barplot"))
                     p1 <- ggplot(gsa_results_total, aes(x=Description, y=count, fill=p.adjust)) +
                       geom_bar(stat="identity") +
                       coord_flip()
                   }
                   
                   incProgress(3/3, detail = "success")
                 })
    p1
  })
  
  
  
}







# Run the application 
shinyApp(ui = ui, server = server)