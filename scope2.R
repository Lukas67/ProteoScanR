
# workflow as per 
# https://uclouvain-cbio.github.io/SCP.replication/articles/SCoPE2.html
# and https://bioconductor.org/packages/release/bioc/vignettes/scp/inst/doc/scp.html

# install scp package

#BiocManager::install("scp")
#BiocManager::install("scpdata")

# read in SCP data

library("scp")
library("ggplot2")
library("magrittr")
library("dplyr")


# read in MS result table
mqScpData <- read.delim("/home/lukas/Desktop/MS-Data/Lukas/mq-run_150223/combined/txt/evidence.txt")

# create annotation file
# this varies upon experimental design
#quantCols <- grep("Reporter.intensity.\\d", colnames(mqScpData), value = T)


# create sample file
# this needs to be done by the researcher
# sampleAnnotation <- as.data.frame(quantCols)
# colnames(sampleAnnotation) <-c("Channel")
# sampleAnnotation$Raw.file = unique(mqScpData$Raw.file)
# 
# 
# samplesA <- c(replicate(6, "Monocytes_A"))
# samplesB <- c(replicate(6, "Monocytes_B"))
# samples <- c(samplesA, samplesB)
# 
# 
# sampleAnnotation$SampleType <- samples

sampleAnnotation = read.delim("/home/lukas/Desktop/MS-Data/Lukas/mq-run_150223/combined/txt/sampleAnnotation_two_groups.txt")


# create QFeature object
scp <- readSCP(featureData = mqScpData,
               colData = sampleAnnotation,
               channelCol = "Channel",
               batchCol = "Raw.file",
               suffix = paste0("_TMT", 1:12),
               removeEmptyCols = TRUE)


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
#nPSMs <- dims(scp)[1, 1]
#nPeps <- dims(scp)[1, 2]

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
                  i = 1:length(rowDataNames(scp)),
                  PEP = "PEP", # by reference the dart_PEP value is used
                  rowDataName = "qvalue_PSMs")

scp <- pep2qvalue(scp,
                  i = 1:length(rowDataNames(scp)),
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

scp <- aggregateFeaturesOverAssays(scp,
                                  i = 1:length(rowDataNames(scp)),
                                  fcol = "Modified.sequence",
                                  name = paste0("peptides_", names(scp)),
                                  fun = matrixStats::colMedians, na.rm = TRUE)
# filter per sample type
# scp <- scp[, scp$SampleType %in% c("Blank", "Macrophage", "Monocyte"), ]

#scp <- joinAssays(scp,
#                  i = 1:length(rowDataNames(scp)),
#                  name = "Peptides")

# filter by median intensity
# calculation of median not applicable if only one run observed

file_name <- sampleAnnotation$Raw.file[1]
peptides <- paste("peptides_", as.character(file_name), sep = "")

medians <- colMedians(assay(scp[[peptides]]), na.rm = TRUE)
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

scp <- medianCVperCell(scp,
                       i = length(rowDataNames(scp)),
                       groupBy = "Leading.razor.protein",
                       nobs = 5, 
                       norm = "div.median",
                       na.rm = TRUE,
                       colDataName = "MedianCV")

getWithColData(scp, peptides) %>%
  colData %>%
  data.frame %>%
  ggplot(aes(x = MedianCV,
             fill = SampleType)) +
  geom_boxplot() +
  geom_vline(xintercept = 0.65)

# filter out all samples with coefficient of variation larger than 
# variable can be adjusted!

scp <- scp[, !is.na(scp$MedianCV) & scp$MedianCV < 0.65, ]

# Normalization of peptide data
# samples --> divide relative intensities by median relative intensities
# peptides --> divide relative intensities by the mean relative intensities

# divide column by median
scp <-normalizeSCP(scp, 2, name="peptides_norm_col", method = "div.median")
#divide rows my mean
scp <-normalizeSCP(scp, 3, name="peptides_norm", method = "div.mean")

# remove peptides with high missing rate

scp <- filterNA(scp,
                i = "peptides_norm",
                pNA = 0.99)
# log-transform 
scp <- logTransform(scp,
                    base = 2,
                    i = "peptides_norm",
                    name = "peptides_log")

# aggregate peptide to protein

scp <- aggregateFeatures(scp,
                         i = "peptides_log",
                         name = "proteins",
                         fcol = "Leading.razor.protein",
                         fun = matrixStats::colMedians, na.rm = TRUE)


# normalization of protein data

## Center columns with median
scp <- sweep(scp, i = "proteins",
             MARGIN = 2,
             FUN = "-",
             STATS = colMedians(assay(scp[["proteins"]]),
                                na.rm = TRUE),
             name = "proteins_norm_col")

## Center rows with mean
scp <- sweep(scp, i = "proteins_norm_col",
             MARGIN = 1,
             FUN = "-",
             STATS = rowMeans(assay(scp[["proteins_norm_col"]]),
                              na.rm = TRUE),
             name = "proteins_norm")

# missing value imputation

# show missing values
scp[["proteins_norm"]] %>%
  assay %>%
  is.na %>%
  mean

# longFormat(scp[, , "proteins_norm"]) %>%
#   data.frame %>%
#   group_by(colname) %>%
#   summarize(missingness = mean(is.na(value))) %>%
#   ggplot(aes(x = missingness)) +
#   geom_histogram()
# 

# use knn to impute missing values

scp <- impute(scp,
              i = "proteins_norm",
              name = "proteins_imptd",
              method = "knn",
              k = 3, rowmax = 1, colmax= 1,
              maxp = Inf, rng.seed = 1234)

# show missing values again
scp[["proteins_imptd"]] %>%
  assay %>%
  is.na %>%
  mean

# batch correction
# upon multiple runs

 sce <- getWithColData(scp, "proteins_norm")

 batch <- colData(sce)$Raw.file
 model <- model.matrix(~ SampleType, data = colData(sce))
 # library(sva)
 # assay(sce) <- ComBat(dat = assay(sce),
 #                      batch = batch,
 #                      mod = model)
 # 

 scp <- addAssay(scp,
                 y = sce,
                 name = "proteins_final")

 scp <- addAssayLinkOneToOne(scp,
                             from = "proteins_norm",
                             to = "proteins_final")



# dimensionality reduction
# changed from reference due no missing values

library(scater)

scp[["proteins_final"]] <- runPCA(scp[["proteins_final"]],
                                   ncomponents = 5,
                                   ntop = Inf,
                                   scale = TRUE,
                                   exprs_values = 1,
                                   name = "PCA")


plotReducedDim(scp[["proteins_final"]],
               dimred = "PCA",
               colour_by = "SampleType",
               point_alpha = 1)

# UMAP
scp[["proteins_final"]] <- runUMAP(scp[["proteins_final"]],
                                   ncomponents = 2,
                                   ntop = Inf,
                                   scale = TRUE,
                                   exprs_values = 1,
                                   n_neighbors = 3,
                                   dimred = "PCA",
                                   name = "UMAP")

plotReducedDim(scp[["proteins_final"]],
               dimred = "UMAP",
               colour_by = "SampleType",
               point_alpha = 1)




## Get the features
subsetByFeature(scp, "Q9ULV4") %>%
  ## Format the `QFeatures` to a long format table
  longFormat(colvars = c("Raw.file", "SampleType", "Channel")) %>%
  data.frame %>%
  ## This is used to preserve ordering of the samples and assays in ggplot2
  mutate(assay = factor(assay, levels = names(scp)),
         Channel = sub("Reporter.intensity.", "Label", Channel),
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

## limma analysis 
# differential expression between two groups

library(limma)

# create expression matrix from protein data
# Log-transformed expression data in a matrix:
#  Each column represents an experiement, and each row represents a detected gene/probe.
# not sure if the data should be aggregated

exp_matrix <- data.frame(assay(scp[["proteins_final"]]))
colnames(exp_matrix) <- scp$SampleType

# exp_matrix <- aggregate(. ~ SampleType, data = data_to_aggregate, FUN = mean, na.rm = TRUE)
# exp_matrix <- data.frame(t(exp_matrix))
# names(exp_matrix) <- exp_matrix[1,]
# exp_matrix <- exp_matrix[-1,]

# create design matrix 
# return the count of individual colnames (Types of experiment)
# Sample data
factor_var <- factor(colnames(exp_matrix))

# Create lookup table
levels <- unique(factor_var)
num_values <- seq_along(levels)
lookup_table <- data.frame(factor_var = levels, num_var = num_values)

# Convert factor variable to numeric
num_var <- sapply(factor_var, function(x) {
  lookup_table$num_var[lookup_table$factor_var == x]
})

# Replace original variable with numeric version
factor_var <- factor(num_var)

# define design without intercept
design <- model.matrix(~0+factor_var)
# assign the column names
colnames(design) <- c(unique(scp$SampleType))

# define design with intercept
design <- model.matrix(~factor_var)
# assign the column names
colnames(design) <- c("Intercept", unique(scp$SampleType)[-1])


# Fit the expression matrix to a linear model
fit <- lmFit(exp_matrix, design)
# Compute contrast
# fit_contrast <- contrasts.fit(fit, cont_matrix)
# Bayes statistics of differential expression
# *There are several options to tweak!*
# fit_contrast <- eBayes(fit_contrast)

fit <- eBayes(fit)
# Generate a vocalno plot to visualize differential expression
#volcanoplot(fit_contrast)

volcanoplot(fit)
# Generate a list of top 100 differentially expressed genes
# top_genes <- topTable(fit_contrast, number = 100, adjust = "BH")
top_genes <- topTable(fit, number = 100, adjust = "fdr")

# Summary of results (number of differentially expressed genes)
# result <- decideTests(fit_contrast)

result <- decideTests(fit)

summary(result)


# QQ-plots for differential expression

ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma

par(mfrow=c(1,2))
qqt(ordinary.t, df=fit$df.residual, main="Ordinary t")
abline(0,1)
qqt(fit$t, df=fit$df.total,main="Moderated t")
abline(0,1)
par(mfrow=c(1,1))


# similar to eBayes but with fold change threshold
treat(fit, lfc=0, trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))



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
