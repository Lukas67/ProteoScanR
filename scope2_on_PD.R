
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
evidence_data <- read.delim("/home/lukas/Desktop/MS-Data/Lukas/Monocytes/gerhard_pd/20230202_HF2_08_Ujjwal_Monocytes_TMT12_PSMs.txt")

meta_data_0 = read.delim("/home/lukas/Desktop/MS-Data/Lukas/Monocytes/gerhard_pd/sampleAnnotate_tabdel.txt")

# define cutoff values for example

PIF_cutoff <- 0.9
qvalue_cutoff <- 0.05
nObs_pep_razrpr <- 5
MedCV_thresh <- 0.65
pNA <- 0.99

if ("Spectrum.File" %in% colnames(evidence_data) && ("Spectrum.File" %in% colnames(meta_data_0))) {
  
  # parse colnames to maxquant style
  colnames(evidence_data)[colnames(evidence_data) == "Spectrum.File"] <- "Raw.file"
  colnames(meta_data_0)[colnames(meta_data_0) == "Spectrum.File"] <- "Raw.file"
  
  if ("Percolator.PEP" %in% colnames(evidence_data)) {
    colnames(evidence_data)[colnames(evidence_data) == "Percolator.PEP"] <- "PEP"
  }
  if ("Master.Protein.Accessions" %in% colnames(evidence_data)) {
    colnames(evidence_data)[colnames(evidence_data) == "Master.Protein.Accessions"] <- "Leading.razor.protein"
  }
  if ("Annotated.Sequence" %in% colnames(evidence_data)) {
    colnames(evidence_data)[colnames(evidence_data) == "Annotated.Sequence"] <- "Modified.sequence"
  }
  
  # calculate PIF for the MS-spectra
  evidence_data$PIF <- evidence_data$Precursor.Intensity / evidence_data$Intensity
  # could also be this
  #evidence_data$PIF <- evidence_data$Master.Scans / evidence_data$Last.Scan
  
  PD_validator <- TRUE  
}


scp_0 <- readSCP(featureData = evidence_data,
                 colData = meta_data_0,
                 channelCol = "Channel",
                 batchCol = "Raw.file",
                 removeEmptyCols = TRUE)

#manual selection of SampleTypes to exclude from the analysis
#scp_0 <- scp_0[, !(scp_0$SampleType %in%  input$selectedSampleType_to_exclude)]      


# clean missing data
scp_0 <- zeroIsNA(scp_0, i=1:length(rowDataNames(scp_0)))


if (PD_validator) {
  scp_0 <- filterFeatures(scp_0,
                          ~ Confidence == "High" &
                            !is.na(PIF) & PIF > PIF_cutoff)
  
} else {
  scp_0 <- filterFeatures(scp_0,
                          ~ Reverse != "+" &
                            Potential.contaminant != "+" &
                            !is.na(PIF) & PIF > PIF_cutoff)
}

scp_0 <- pep2qvalue(scp_0,
                    i = names(scp_0),
                    PEP = "PEP", # by reference the dart_PEP value is used
                    rowDataName = "qvalue_PSMs")

scp_0 <- pep2qvalue(scp_0,
                    i = names(scp_0),
                    PEP = "PEP",
                    groupBy = "Leading.razor.protein",
                    rowDataName = "qvalue_proteins")

scp_0 <- filterFeatures(scp_0, ~ qvalue_proteins < qvalue_cutoff)

scp_0 <- aggregateFeaturesOverAssays(scp_0,
                                     i = names(scp_0),
                                     fcol = "Modified.sequence",
                                     name = paste0("peptides_", names(scp_0)),
                                     fun = matrixStats::colMedians, na.rm = TRUE)

if (length(unique(meta_data_0$Raw.file)) > 1) {
  scp_0 <- joinAssays(scp_0,
                      i = ((length(names(scp_0))/2)+1):length(names(scp_0)),
                      name = "peptides")
}

# calculate median reporter IO intensity
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

scp_0 <- scp_0[, !is.na(scp_0$MedianCV) & scp_0$MedianCV < MedCV_thresh, ]

if (length(peptide_file) > 1) {
  scp_0 <- filterNA(scp_0,
                    i = "peptides",
                    pNA = pNA)
} else {
  scp_0 <- filterNA(scp_0,
                    i = peptide_file,
                    pNA = pNA)
}

if (length(peptide_file) > 1) {
  if (PD_validator) {
    
  } else {
    scp_0 <- aggregateFeatures(scp_0,
                               i = "peptides",
                               name = "proteins",
                               fcol = "Leading.razor.protein",
                               fun = matrixStats::colMedians, na.rm = TRUE)    
  }

} else {
  if (PD_validator) {
    sce <- getWithColData(scp_0, peptide_file)
    rownames(sce) <- NULL 
    
    scp_0 <- aggregateFeatures(sce,
                               i = peptide_file,
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







