
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
akos_data <- read.csv("/home/lukas/Downloads/20230202_HF2_08_Ujjwal_Monocytes_TMT12_proteins.csv")



# create annotation file
# this varies upon experimental design
quantCols <- grep("Reporter.intensity.\\d", colnames(mqScpData), value = T)


# create sample file
# this needs to be done by the researcher
sampleAnnotation <- as.data.frame(quantCols)
colnames(sampleAnnotation) <-c("Channel")
sampleAnnotation$Raw.file = unique(mqScpData$Raw.file)


samplesA <- c(replicate(6, "Monocytes-A"))
samplesB <- c(replicate(6, "Monocytes-B"))
samples <- c(samplesA, samplesB)


sampleAnnotation$SampleType <- samples

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
                        !is.na(PIF) & PIF > 0.8)

# filter runs with too few features
# however only one run is done
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

library(stringr)

find_peptide_list <- function(scp_object) {
  scp_cols <- colnames(scp_object)
  peptide_list <- str_extract(scp_cols, "^peptides_.+$")
  return(peptide_list)
}

peptide_list <- find_peptide_list(scp)

medians <- colMedians(assay(scp[[]]), na.rm = TRUE)
scp$MedianRI <- medians


colData(scp) %>%
  data.frame %>%
  ggplot() +
  aes(x = MedianRI, 
      y = SampleType,
      fill = SampleType) +
  geom_boxplot() +
  scale_x_log10()