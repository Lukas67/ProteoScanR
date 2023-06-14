
# workflow as per 
# https://uclouvain-cbio.github.io/SCP.replication/articles/SCoPE2.html
# and https://bioconductor.org/packages/release/bioc/vignettes/scp/inst/doc/scp.html

# install scp package

#BiocManager::install("scp")
#BiocManager::install("scpdata")

# read in SCP data

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
library("tibble")

# read in MS result table
mqScpData <- read.delim("/home/lukas/Desktop/MS-Data/Lukas/Jun13/combined/txt/evidence.txt")

sampleAnnotation = read.delim("/home/lukas/Desktop/MS-Data/Lukas/Jun13/combined/txt/sampleAnnotation.txt")



 
# # create sample annotation
# Raw.files <- unique(mqScpData$Raw.file)
# 
# library(stringr)
# Raw.files <- Raw.files[str_order(str_extract_all(Raw.files, "\\d+"), numeric = TRUE)]
# sampleAnnotation <- data.frame(Raw.file=rep(Raw.files, each=12))
# 
# channels <- colnames(mqScpData[ , grepl( "Reporter.intensity.corrected" , names(mqScpData) ) ])
# sampleAnnotation$Channel <- rep(channels, 12)
# 
# sampletypes <- c(rep("HC", 4), rep("HIV_MetS", 4), rep("HIV_noMetS", 4))
# sampleAnnotation$SampleType <- rep(sampletypes, 12)
# 
# outfile <- "/home/lukas/Desktop/MS-Data/Lukas/Jun13/combined/txt/sampleAnnotation.txt"
# 
# library("readr")
# write_delim(sampleAnnotation, outfile, delim = "\t")


# create QFeature object
scp <- readSCP(featureData = mqScpData,
               colData = sampleAnnotation,
               channelCol = "Channel",
               batchCol = "Raw.file",
               removeEmptyCols = TRUE)



# skip pool sample
if ("Pool" %in% sampleAnnotation$SampleType) {
  scp <- scp[, scp$SampleType !=  "Pool"]
}



# plot runs (not needed - sample unfractionized)
# only one run
# plot(scp)

# clean missing data
scp <- zeroIsNA(scp, i=1:length(rowDataNames(scp)))

# filter PSM
# filter out potential contaminants
# filter out matches to decoy database
# keep PSMs with high PIF (parental ion fraction)
scp <- filterFeatures(scp,
                      ~ Reverse != "+" &
                        Potential.contaminant != "+" &
                        !is.na(PIF) & PIF > 0.10)


nPSMs <- dims(scp)[1, ]

ggplot(data.frame(nPSMs)) +
  aes(x = nPSMs) +
  geom_histogram(binwidth = 50) +
  geom_vline(xintercept = 500)


# filter runs with too few features
# however only one run is done --> changes nothing
keepAssay <- dims(scp)[1, ] > 150
scp <- scp[, , keepAssay]

# filter out based on SCP metric
# sample to carrier ratio --> RI of single cell sample / RI carrier channel
# carrier channel is meant to boost peptide id. rate


#however no carrier channel is used
#scp <- computeSCR(scp,
#                  i = 1:3,
#                  colvar = "SampleType",
#                  carrierPattern = "Carrier",
#                  samplePattern = "Macrophage|Monocyte",
#                  sampleFUN = "mean",
#                  rowDataName = "MeanSCR")

#rbindRowData(scp, i = 1:3) %>%
#  data.frame %>%
#  ggplot(aes(x = MeanSCR)) +
#  geom_histogram() +
#  geom_vline(xintercept = c(1/200, 0.1),
#             lty = c(2, 1)) +
#  scale_x_log10()

#scp <- filterFeatures(scp,
#                      ~ !is.na(MeanSCR) &
#                        MeanSCR < 0.1)
#> 'MeanSCR' found in 3 out of 3 assay(s)





# compute qvalues_PSMs to filter out by FDR

scp <- pep2qvalue(scp,
                  i = names(scp),
                  PEP = "PEP", # by reference the dart_PEP value is used
                  rowDataName = "qvalue_PSMs")

scp <- pep2qvalue(scp,
                  i = names(scp),
                  PEP = "PEP",
                  groupBy = "Leading.razor.protein",
                  rowDataName = "qvalue_proteins")


scp <- filterFeatures(scp, ~ qvalue_proteins < 0.01)



# correct in between run variation 
# when only one run observed --> not needed

#scp <- divideByReference(scp,
#                         i = 1:length(rowDataNames(scp)),
#                         colvar = "SampleType",
#                         samplePattern = ".",
#                         refPattern = "Reference")

# aggregate PSMs to peptides

# can be called in the function argument of aggregateFeaturesOverAssays
remove.duplicates <- function(x)
  apply(x, 2, function(xx) xx[which(!is.na(xx))[1]] )

# matrixStats::colMedians, na.rm = TRUE aggregate over median 
scp <- aggregateFeaturesOverAssays(scp,
                                  i = names(scp),
                                  fcol = "Modified.sequence",
                                  name = paste0("peptides_", names(scp)),
                                  fun = matrixStats::colMedians, na.rm = TRUE)


# join assays to one
if (length(names(scp)) > 1) {
  scp <- joinAssays(scp,
                    i = ((length(names(scp))/2)+1):length(names(scp)),
                    name = "peptides")
}
  

# filter by median intensity
file_name <- unique(sampleAnnotation$Raw.file)
peptide_file <- paste("peptides_", as.character(file_name), sep = "")

if (length(peptide_file) > 1) {
  medians <- c()
  for (assay_name in peptide_file) {
    new_medians <- colMedians(assay(scp[[assay_name]]), na.rm = TRUE)
    medians <- c(medians, new_medians)
  }
} else {
  medians <- colMedians(assay(scp[[peptide_file]]), na.rm = TRUE)
}
scp$MedianRI <- medians

colData(scp) %>%
  data.frame %>%
  ggplot() +
  aes(x = MedianRI, 
      y = SampleType,
      fill = SampleType) +
  geom_boxplot() +
  scale_x_log10()

# calculate median CV (coefficient of variation)
if (length(peptide_file) > 1) {
  scp <- medianCVperCell(scp,
                         i = "peptides",
                         groupBy = "Leading.razor.protein",
                         nobs = 5, 
                         norm = "div.median",
                         na.rm = TRUE,
                         colDataName = "MedianCV")
} else {
  scp <- medianCVperCell(scp,
                         i = peptide_file,
                         groupBy = "Leading.razor.protein",
                         nobs = 5, 
                         norm = "div.median",
                         na.rm = TRUE,
                         colDataName = "MedianCV")
}

if (length(peptide_file) > 1) {
  getWithColData(scp, "peptides") %>%
    colData %>%
    data.frame %>%
    ggplot(aes(x = MedianCV,
               fill = SampleType)) +
    geom_boxplot() +
    geom_vline(xintercept = 0.65)  
} else {
  getWithColData(scp, peptide_file) %>%
    colData %>%
    data.frame %>%
    ggplot(aes(x = MedianCV,
               fill = SampleType)) +
    geom_boxplot()
}


# filter out all samples with coefficient of variation larger than 
# variable can be adjusted!

scp <- scp[, !is.na(scp$MedianCV) & scp$MedianCV < 0.65, ]

# Normalization of peptide data
# samples --> divide relative intensities by median relative intensities
# peptides --> divide relative intensities by the mean relative intensities

# divide column by median
# scp <-normalizeSCP(scp, 2, name="peptides_norm_col", method = "div.median")
# #divide rows my mean
# scp <-normalizeSCP(scp, 3, name="peptides_norm", method = "div.mean")

# remove peptides with high missing rate



if (length(peptide_file) > 1) {
  scp <- filterNA(scp,
                  i = "peptides",
                  pNA = 0.99)
} else {
  scp <- filterNA(scp,
                  i = peptide_file,
                  pNA = 0.99)
}


# # sqrt transformation
# scp <- sweep(scp, i="peptides_norm",
#              MARGIN = 2,
#              FUN="^",
#              STATS=2,
#              name="peptide_sqrt")
# 
# #log-transform
# scp <- logTransform(scp,
#                     base = 2,
#                     i = "peptides_norm",
#                     name = "peptides_log")
# 

# aggregate peptide to protein
if (length(peptide_file) > 1) {
  scp <- aggregateFeatures(scp,
                           i = "peptides",
                           name = "proteins",
                           fcol = "Leading.razor.protein",
                           fun = matrixStats::colMedians, na.rm = TRUE)
  
} else {
  scp <- aggregateFeatures(scp,
                           i = peptide_file,
                           name = "proteins",
                           fcol = "Leading.razor.protein",
                           fun = matrixStats::colMedians, na.rm = TRUE)
}


#log-transform
# scp <- logTransform(scp,
#                     base = 2,
#                     i = "proteins",
#                     name = "proteins_transf")

# box cox transf
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


protein_matrix <- assay(scp[["proteins"]])

prot_lm <- lm(protein_matrix ~ 1)

b <- boxcox_1(prot_lm)

b




# Exact lambda
lambda <- b$x[which.max(b$y)]
lambda

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


sce <- getWithColData(scp, "proteins")

scp <- addAssay(scp,
                y = sce,
                name = "proteins_transf")

scp <- addAssayLinkOneToOne(scp,
                            from = "proteins",
                            to = "proteins_transf")

assay(scp[["proteins_transf"]]) <- protein_matrix

# normalization of protein data

## SCoPE2 normalization

# # Center columns with median
#  scp <- sweep(scp, i = "proteins",
#               MARGIN = 2,
#               FUN = "-",
#               STATS = colMedians(assay(scp[["proteins"]]),
#                                  na.rm = TRUE),
#               name = "proteins_norm_col")
# 
# ## Center rows with mean
#  scp <- sweep(scp, i = "proteins_norm_col",
#               MARGIN = 1,
#               FUN = "-",
#               STATS = rowMeans(assay(scp[["proteins_norm_col"]]),
#                                na.rm = TRUE),
#               name = "proteins_norm")

 
 # CONSTANd normalization method

# # make MA plot first
# MAplot <- function(x,y,use.order=FALSE, R=NULL, cex=1.6, showavg=TRUE) {
#   # catch unequal size of matrices
#   if (dim(x)[2] != dim(y)[2]) {
#     if (dim(x)[2] > dim(y)[2]) {
#       x <- x[, 1:dim(y)[2]]
#     } else if (dim(y)[2] > dim(x)[2]) {
#       y <- y[, 1:dim(x)[2]]
#     }
#   }
#   
#   # make an MA plot of y vs. x that shows the rolling average,
#   M <- log2(y/x)
#   xlab = 'A'
#   if (!is.null(R)) {r <- R; xlab = "A (re-scaled)"} else r <- 1
#   A <- (log2(y/r)+log2(x/r))/2
#   if (use.order) {
#     orig.order <- order(A)
#     A <- orig.order
#     M <- M[orig.order]
#     xlab = "original rank of feature magnitude within IPS"
#   }
#   # select only finite values
#   use <- is.finite(M)
#   A <- A[use]
#   M <- M[use]
#   # plot
#   print(var(M))
#   plot(A, M, xlab=xlab, cex.lab=cex, cex.axis=cex)
#   # rolling average
#   if (showavg) { lines(lowess(M~A), col='red', lwd=5) }
# }
# 
# MAplot(assay(scp[["proteins_transf"]])[,1:6], assay(scp[["proteins_transf"]])[,7:12])
# 
# #find all pairwise indeces
# st_indeces <- split(seq_along(scp$SampleType), scp$SampleType)
# 
# comp_list <- apply(combn(unique(scp$SampleType),2),2,paste, collapse="-")
# 
# index_combis <- apply(combn(st_indeces,2),2,paste)
# 
# user_choice <- comp_list[1]
# user_choice_vector <- strsplit(user_choice, split = "-")
# choice_A <- user_choice_vector[[1]][1]
# choice_B <- user_choice_vector[[1]][2]
# index_A <- st_indeces[choice_A]
# index_B <- st_indeces[choice_B]
# 
# MAplot(assay(scp[["proteins_transf"]][,index_A[[1]]]), assay(scp[["proteins_transf"]][,index_B[[1]]]))
# 
# protein_matrix <- assay(scp[["proteins_transf"]])
# protein_matrix <- CONSTANd(protein_matrix)
# 
# sce <- getWithColData(scp, "proteins")
# 
# scp <- addAssay(scp,
#                 y = sce,
#                 name = "proteins_norm")
# 
# scp <- addAssayLinkOneToOne(scp,
#                             from = "proteins_transf",
#                             to = "proteins_norm")
# 
# #assay(scp[["proteins_norm"]]) <- assay(scp[["proteins_transf"]])
# assay(scp[["proteins_norm"]]) <- protein_matrix$normalized_data

protein_matrix <- assay(scp[["proteins_transf"]])

sce <- getWithColData(scp, "proteins_transf")

scp <- addAssay(scp,
                  y = sce,
                  name = "proteins_norm")

scp <- addAssayLinkOneToOne(scp,
                              from = "proteins_transf",
                              to = "proteins_norm")

protein_matrix <- normalizeQuantiles(protein_matrix)

assay(scp[["proteins_norm"]]) <- protein_matrix




# missing value imputation
# show missing values
scp[["proteins_norm"]] %>%
  assay %>%
  is.na %>%
  mean

longFormat(scp[, , "proteins_norm"]) %>%
  data.frame %>%
  group_by(colname) %>%
  summarize(missingness = mean(is.na(value))) %>%
  ggplot(aes(x = missingness)) +
  geom_histogram()


# use knn to impute missing values
library(impute)
library(sva)

if (length(peptide_file) > 1) {
  # scp <- impute(scp,
  #               i = "proteins_norm",
  #               name = "proteins_imptd",
  #               method = "knn",
  #               k = 3, rowmax = 1, colmax= 1,
  #               maxp = Inf, rng.seed = as.numeric(gsub('[^0-9]', '', Sys.Date())))
  protein_matrix <- assay(scp[["proteins_norm"]])
  
  sce <- getWithColData(scp, "proteins_norm")
  
  scp <- addAssay(scp,
                  y = sce,
                  name = "proteins_imptd")
  
  scp <- addAssayLinkOneToOne(scp,
                              from = "proteins_norm",
                              to = "proteins_imptd")
  
  protein_matrix <- impute.knn(protein_matrix, 
                               k=3, 
                               rowmax = 1, 
                               colmax = 1, 
                               maxp = Inf, 
                               rng.seed = as.numeric(gsub('[^0-9]', '', Sys.Date())))
  
  
  assay(scp[["proteins_imptd"]]) <- protein_matrix$data
  } 


  
  sce <- getWithColData(scp, "proteins_imptd")
  
  batch <- colData(sce)$Raw.file
  model <- model.matrix(~SampleType, data = colData(sce))
  
  Combat_batchC <- function(i_exp_matrix, i_batch, i_model) {
    out <- tryCatch(
      {
        i_exp_matrix <- ComBat(dat = i_exp_matrix,
                               batch = i_batch,
                               mod = i_model)
        print("batch corrected and optimized")
      } ,
      error = function(cond) {
        i_exp_matrix <- ComBat(dat = i_exp_matrix,
                               batch = i_batch)
        print("confounder detected! just batch corrected")
      }
    )
    return(i_exp_matrix)
  }
  
  assay(sce) <- Combat_batchC(assay(sce), batch, model) 
  
  
  
  scp <- addAssay(scp,
                    y = sce,
                    name = "proteins_dim_red")
  
  scp <- addAssayLinkOneToOne(scp,
                                from = "proteins_imptd",
                                to = "proteins_dim_red")
  

# show missing values again
scp[["proteins_imptd"]] %>%
  assay %>%
  is.na %>%
  mean

# batch correction
# upon multiple runs


# dimensionality reduction
# changed from reference due no missing values

library(scater)

scp[["proteins_dim_red"]] <- runPCA(scp[["proteins_dim_red"]],
                                   ncomponents = 5,
                                   ntop = Inf,
                                   scale = TRUE,
                                   exprs_values = 1,
                                   name = "PCA")


plotReducedDim(scp[["proteins_dim_red"]],
               dimred = "PCA",
               colour_by = "Channel",
               point_alpha = 1,
               )

# UMAP
scp[["proteins_dim_red"]] <- runUMAP(scp[["proteins_dim_red"]],
                                   ncomponents = 3,
                                   ntop = Inf,
                                   scale = TRUE,
                                   exprs_values = 1,
                                   n_neighbors = 3,
                                   dimred = "PCA",
                                   name = "UMAP")

plotReducedDim(scp[["proteins_dim_red"]],
               dimred = "UMAP",
               colour_by = "SampleType",
               point_alpha = 1)


channelstring <- gsub("[0-9]{1,2}$","", scp$Channel[1]) 


## Get the features
subsetByFeature(scp, "E9PAV3") %>%
  ## Format the `QFeatures` to a long format table
  longFormat(colvars = c("Raw.file", "SampleType", "Channel")) %>%
  data.frame %>%
  ## This is used to preserve ordering of the samples and assays in ggplot2
  mutate(assay = factor(assay, levels = names(scp)),
         Channel = sub(channelstring, "", Channel),
  ) %>%
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



## limma analysis 
# differential expression between two groups

library(limma)
# create expression matrix from protein data
# Log-transformed expression data in a matrix:
#  Each column represents an experiement, and each row represents a detected gene/probe.


# create expression matrix
exp_matrix <- data.frame(assay(scp[["proteins_dim_red"]]))
colnames(exp_matrix) <- scp$SampleType

# check if the data is normal distributed
sample_types <- unique(colnames(exp_matrix))

qqnorm(exp_matrix[[sample_types[3]]], main = paste("QQ-plot of ", sample_types[3]))
qqline(exp_matrix[[sample_types[3]]])

hist(exp_matrix[[sample_types[3]]])
# create design matrix 
# return the count of individual colnames (Types of experiment)
# Sample data
# factor_var <- factor(colnames(exp_matrix))
# 
# factorize_var <- function(i_vector) {
#   fact_vect <- factor(i_vector)
#   levels <- unique(fact_vect)
#   num_values <- seq_along(levels)
#   lookup_table <- data.frame(fact_vect = levels, num_var = num_values)
#   
#   num_var <- sapply(fact_vect, function(x) {
#     lookup_table$num_var[lookup_table$fact_vect == x]
#   })
#   
#   fact_vect <- factor(num_var)
#   return(fact_vect)
# }
# 
# factors_sample_type <- factorize_var(colnames(exp_matrix))
# 
# user_input <- c(1,2,3,4,5,6,1,2,3,4,5,6)
# user_chosen_name <- "patient"
# 
# factors_patient <- factorize_var(user_input)
# 
# patient <- c(1,2,3,4,5,6,1,2,3,4,5,6)
# 
# design_frame <- data.frame(cbind(patient, scp$SampleType))
# # define design without intercept
# design <- model.matrix(~ 0 +  ., data=design_frame)
# # assign the column names
# user_colnames <- sprintf(paste(as.character(user_chosen_name),"[%s]"), seq(2:length(unique(factors_patient)))+1)
# 
# colnames(design) <- c(unique(scp$SampleType), user_colnames)



# # Fit the expression matrix to a linear model
# fit <- lmFit(exp_matrix, design)
# # Compute contrast
# # fit_contrast <- contrasts.fit(fit, cont_matrix)
# # Bayes statistics of differential expression
# # *There are several options to tweak!*
# # fit_contrast <- eBayes(fit_contrast)
# 
# fit <- eBayes(fit)
# # Generate a vocalno plot to visualize differential expression
# #volcanoplot(fit_contrast)
# 
# volcanoplot(fit)
# 
# 
# # Generate a list of top 100 differentially expressed genes
# # top_genes <- topTable(fit_contrast, number = 100, adjust = "BH")
# top_genes <- topTable(fit, number = 100, adjust = "fdr")



# Summary of results (number of differentially expressed genes)
# result <- decideTests(fit_contrast)

# result <- decideTests(fit)
# 
# summary(result)


# QQ-plots for differential expression

# ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
# 
# par(mfrow=c(1,2))
# qqt(ordinary.t, df=fit$df.residual, main="Ordinary t")
# abline(0,1)
# qqt(fit$t, df=fit$df.total,main="Moderated t")
# abline(0,1)
# par(mfrow=c(1,1))
# 
# 
# # similar to eBayes but with fold change threshold
# treat(fit, lfc=0, trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))



# define design matrix calculation for 4 different effects

# --> Additive effect
# --> Interaction effect
# --> Nested Factor
# --> mixed effect models


## Studies with one comparative kind of explanatory variables
# factor var need to be changed to a covariate as well

# use intercept when comparing a treatment versus control
# first column will be taken into intercept variable (0 point)
# design <- model.matrix(~factor_var)

# when performing pairwise comparisons
# group several variables to one factor --> simplifies experiment
# design <- model.matrix(~0+factor_var)

## two or multiple treatment variables can be modelled with different approaches

# additive model
# design <- model.matrix(~factor_var + treatment)

# interactive model
# design <- model.matrix(~factor_var * treatment)

## time models 
# combination of time and treatment
# design <- model.matrix(~factor_var * time)

# linear series
# design <- model.matrix(~time)




# defining the linear model 

# fit <- lmFit(exp_matrix ~ design)

# in order to define 0 as the y intercept redefine design matrix


# read the object
scp_0 <- scp

# fetch user input for contrasts
# selectedComp_stat <- c("A", "B", "D")
# # update dataframe according to user selection
# scp_0 <- scp_0[, scp_0$SampleType %in%  selectedComp_stat]
# update expression matrix according to user selection
exp_matrix_0 <- assay(scp_0[["proteins_dim_red"]])

# # fetch multiple factors for analysis
# user_choice <- c("Raw.file")#, "Technician")
# #user_choice <- c("Technician")
# 
# # fetch the metadata to factorize
# fetched_factor <- colData(scp_0)[user_choice]
# design_frame <- cbind(fetched_factor, scp_0$SampleType)
# design <- model.matrix(~0 + . , data=design_frame)

# single factor analysis
design <- model.matrix(~0 + scp_0$SampleType)
colnames(design)<-unique(scp_0$SampleType)
fit <- lmFit(exp_matrix_0, design)
cont.matrix <- makeContrasts(contrasts = "S_48h-PreOp" , levels=design)
fit <- contrasts.fit(fit, cont.matrix)

fit <- eBayes(fit)

results <- decideTests(fit)
genetable <- topTable(fit, number = Inf, adjust = "BH")

volcanoplot(fit)


library(plotly)

# load data.frame
toptable <- as.data.frame(genetable)
toptable <- cbind(rownames(toptable), data.frame(toptable, row.names=NULL))
colnames(toptable)[1] <- "protein"

# significance cutoffs
FCcutoff <- 0.05
pCutoff <- 0.05

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




heatmap(cor(t(assay(scp_0[["proteins_dim_red"]]))))

dimred <- data.frame(reducedDim(scp_0[["proteins_dim_red"]], "PCA"))
plot_ly(x=dimred$PC1,
        y=dimred$PC2,
        z=dimred$PC3,
        type="scatter3d",
        mode="markers",
        color=scp$Raw.file)



## validate normalization
library("infotheo")

# select data according to the procedure undertaken
data_final <- assay(scp[["proteins_dim_red"]])
data_raw <- assay(scp[["proteins"]])

# select data according to groups to observe mutual information within
data_final <- data_final[, which(scp$SampleType == "PreOp")]
data_raw <- data_raw[, which(scp$SampleType == "PreOp")]
            
# discretize by equal frequencies
data_final <- discretize(data_final)
data_raw <- discretize(data_raw)


# mutual information is returned in nats
mi_final <- mutinformation(data_final)
mi_raw <- mutinformation(data_raw)

gain <- mi_final-mi_raw
boxplot(stack(gain)$value)

mi_final_stack <- data.frame(stack(mi_final))
mi_raw_stack <- data.frame(stack(mi_raw))

mi_total <- cbind(mi_final_stack$value, mi_raw_stack$value)
colnames(mi_total) <- c("MI_final", "MI_raw")

to_plot <- melt(mi_total)

ggplot(to_plot, aes(x=Var2, y=value)) +
  geom_boxplot(aes(fill=Var2)) +
  ggtitle(paste("Change in mutual information after transformation:", 
                "transform_base_val()", 
                "and norm_method:", "norm_method_val()")) +
  xlab(paste("input$selectedSampleType_normval")) +
  ylab("Change in natural unit of information") + 
  guides(fill=guide_legend(title=""))


# inter group comparision of information content before and after the processes 

data_final <- assay(scp[["proteins_dim_red"]])
data_raw <- assay(scp[["proteins"]])

user_choice <- "PreOP-S_24h"
# split user choice of comp back to sample types
user_choice_vector <- strsplit(user_choice, split = "-")
# and assign them to a variable
choice_A <- user_choice_vector[[1]][1]
choice_B <- user_choice_vector[[1]][2]


data_final_group1 <- discretize(data_final[, which(scp$SampleType == choice_A)])
data_raw_group1 <- discretize(data_raw[, which(scp$SampleType == choice_A)])

data_final_group2 <- discretize(data_final[, which(scp$SampleType == choice_B)])
data_raw_group2 <- discretize(data_raw[, which(scp$SampleType == choice_B)])

colnames(data_final_group1) <- paste( colnames(data_final_group1), "grp1")
colnames(data_final_group2) <- paste( colnames(data_final_group2), "grp2")

colnames(data_raw_group1) <- paste( colnames(data_raw_group1), "grp1")
colnames(data_raw_group2) <- paste( colnames(data_raw_group2), "grp2")

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

mi_total <- data.frame(mi_final=mi_final, mi_raw=mi_raw)

to_plot <- melt(mi_total)

ggplot(to_plot, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  ggtitle(paste("Change in mutual information after transformation:", 
                "transform_base_val()", 
                "and norm_method:", "norm_method_val()")) +
  xlab(paste("Comparison between sample types")) +
  ylab("Change in natural unit of information") + 
  guides(fill=guide_legend(title=""))
