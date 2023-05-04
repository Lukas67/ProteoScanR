
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
mqScpData <- read.delim("/home/lukas/Desktop/MS-Data/Lukas/Monocytes/gerhard_pd/20230202_HF2_08_Ujjwal_Monocytes_TMT12_PSMs.txt")

sampleAnnotation = read.delim("/home/lukas/Desktop/MS-Data/Lukas/Monocytes/gerhard_pd/sampleAnnotate_tabdel.txt")

# attempting to calculate PIF
#The calculation involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window
#and is performed before and after the MS/MS scan of interest and interpolated at the recorded time of the MS/MS acquisition

# exact same values? 




# create QFeature object
scp <- readSCP(featureData = mqScpData,
               colData = sampleAnnotation,
               channelCol = "Channel",
               batchCol = "Spectrum.File",
               removeEmptyCols = TRUE)

# clean missing data
scp <- zeroIsNA(scp, i=1:length(rowDataNames(scp)))

scp <- filterFeatures(scp,
                      ~ Confidence == "High")

nPSMs <- dims(scp)[1, ]

ggplot(data.frame(nPSMs)) +
  aes(x = nPSMs) +
  geom_histogram(binwidth = 50) +
  geom_vline(xintercept = 500)

keepAssay <- dims(scp)[1, ] > 150
scp <- scp[, , keepAssay]

scp <- pep2qvalue(scp,
                  i = names(scp),
                  PEP = "Percolator.PEP", # by reference the dart_PEP value is used
                  rowDataName = "qvalue_PSMs")

scp <- pep2qvalue(scp,
                  i = names(scp),
                  PEP = "Percolator.PEP",
                  groupBy = "Master.Protein.Accessions",
                  rowDataName = "qvalue_proteins")


scp <- filterFeatures(scp, ~ qvalue_proteins < 0.01)

scp <- aggregateFeaturesOverAssays(scp,
                                   i = names(scp),
                                   fcol = "Annotated.Sequence",
                                   name = paste0("peptides_", names(scp)),
                                   fun = matrixStats::colMedians, na.rm = TRUE)

if (length(unique(sampleAnnotation$Spectrum.File)) > 1) {
  scp_0 <- joinAssays(scp_0,
                      i = ((length(names(scp))/2)+1):length(names(scp)),
                      name = "peptides")
}

# filter by median intensity
file_name <- unique(sampleAnnotation$Spectrum.File)
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
                         groupBy = "Master.Protein.Accessions",
                         nobs = 5, 
                         norm = "div.median",
                         na.rm = TRUE,
                         colDataName = "MedianCV")
} else {
  scp <- medianCVperCell(scp,
                         i = peptide_file,
                         groupBy = "Master.Protein.Accessions",
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

scp <- scp[, !is.na(scp$MedianCV) & scp$MedianCV < 0.65, ]

if (length(peptide_file) > 1) {
  scp <- filterNA(scp,
                  i = "peptides",
                  pNA = 0.99)
} else {
  scp <- filterNA(scp,
                  i = peptide_file,
                  pNA = 0.99)
}

if (length(peptide_file) > 1) {
  scp <- aggregateFeatures(scp,
                           i = "peptides",
                           name = "proteins",
                           fcol = "Master.Protein.Accessions",
                           fun = matrixStats::colMedians, na.rm = TRUE)
  
} else {
  scp <- aggregateFeatures(scp,
                           i = peptide_file,
                           name = "proteins",
                           fcol = "Master.Protein.Accessions",
                           fun = matrixStats::colMedians, na.rm = TRUE)
}

