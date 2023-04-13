
library("scp")
library("ggplot2")
library("magrittr")
library("dplyr")
library("CONSTANd")
library("MASS")


# read in MS result table
mqScpData <- read.delim("/home/lukas/Desktop/MS-Data/Lukas/mq-run_150223/combined/txt/evidence.txt")


sampleAnnotation = read.delim("/home/lukas/Desktop/MS-Data/Lukas/mq-run_150223/combined/txt/sampleAnnotation_two_groups.txt")

# create QFeature object
scp <- readSCP(featureData = mqScpData,
               colData = sampleAnnotation,
               channelCol = "Channel",
               batchCol = "Raw.file",
               suffix = paste0("_TMT", 1:12),
               removeEmptyCols = TRUE)

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


# Plot number of peptide sequence matches
nPSMs <- dims(scp)[1, ]
#nPSMs <- dims(scp)[1, 1]
#nPeps <- dims(scp)[1, 2]

ggplot(data.frame(nPSMs)) +
  aes(x = nPSMs) +
  geom_histogram(binwidth = 50) +
  geom_vline(xintercept = 500)

# filter runs with lack of features
keepAssay <- dims(scp)[1, ] > 150
scp <- scp[, , keepAssay]


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


# aggregate PSMs to peptides
scp <- aggregateFeaturesOverAssays(scp,
                                   i = 1:length(rowDataNames(scp)),
                                   fcol = "Modified.sequence",
                                   name = paste0("peptides_", names(scp)),
                                   fun = matrixStats::colMedians, na.rm = TRUE)



# filter by median intensity
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

# remove peptides with high missing rate
scp <- filterNA(scp,
                i = peptides, #_norm",
                pNA = 0.99)

# aggregate peptide to protein
scp <- aggregateFeatures(scp,
                         i = peptides,
                         name = "proteins",
                         fcol = "Leading.razor.protein",
                         fun = matrixStats::colMedians, na.rm = TRUE)

#no transform
sce <- getWithColData(scp, "proteins")

scp <- addAssay(scp,
                y = sce,
                name = "proteins_transf")

scp <- addAssayLinkOneToOne(scp,
                            from = "proteins",
                            to = "proteins_transf")

# no normalization applied
sce <- getWithColData(scp, "proteins_transf")

scp <- addAssay(scp,
                y = sce,
                name = "proteins_norm")

scp <- addAssayLinkOneToOne(scp,
                            from = "proteins",
                            to = "proteins_norm")

# missing value imputation
# show missing values
scp[["proteins_transf"]] %>%
  assay %>%
  is.na %>%
  mean



sce <- getWithColData(scp, "proteins_transf")

scp <- addAssay(scp,
                y = sce,
                name = "proteins_final")

scp <- addAssayLinkOneToOne(scp,
                            from = "proteins_transf",
                            to = "proteins_final")


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



library(limma)
# create expression matrix
exp_matrix <- data.frame(assay(scp[["proteins_final"]]))
colnames(exp_matrix) <- scp$SampleType

# check if the data is normal distributed
sample_types <- unique(colnames(exp_matrix))

qqnorm(exp_matrix[[sample_types[1]]], main = paste("QQ-plot of ", sample_types[1]))
qqline(exp_matrix[[sample_types[1]]])

hist(exp_matrix[[sample_types[1]]])


qqnorm(exp_matrix[[sample_types[2]]], main = paste("QQ-plot of ", sample_types[2]))
qqline(exp_matrix[[sample_types[2]]])

hist(exp_matrix[[sample_types[2]]])



# fetch multiple factors for analysis
user_choice <- c("Pair")#, "Technician")
#user_choice <- c("Technician")

# fetch the metadata to factorize
fetched_factor <- colData(scp)[user_choice]
design_frame <- cbind(fetched_factor, scp$SampleType)
design <- model.matrix(~0 + . , data=design_frame)

fit <- lmFit(exp_matrix, design)
fit <- eBayes(fit)

results <- decideTests(fit)
topTable(fit)#, coef = "`scp$SampleType`Monocytes_B")

volcanoplot(fit)


