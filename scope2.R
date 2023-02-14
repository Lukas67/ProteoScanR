# install scp package

#BiocManager::install("scp")
#BiocManager::install("scpdata")

# read in SCP data

library("scp")
library("ggplot2")
library("magrittr")
library("dplyr")


# read in MS result table
mqScpData <- read.delim("/home/lukas/Desktop/MS-Data/Lukas/combined/txt/evidence.txt")


# create annotation file
# this varies upon experimental design
quantCols <- grep("Reporter.intensity.\\d", colnames(mqScpData), value = T)

# create sample file
# this needs to be done by the researcher
sampleAnnotation <- as.data.frame(quantCols)
colnames(sampleAnnotation) <-c("Channel")
sampleAnnotation$Raw.file = unique(mqScpData$Raw.file)


# for example generated --> fix accordingly
generate_var = function(prefix, length) {
  paste(prefix, sprintf('%03d', seq_len(length)), sep = '_')
}

sampleAnnotation$SampleType <- c(generate_var("example", nrow(sampleAnnotation)))


scp <- readSCP(featureData = mqScpData,
               colData = sampleAnnotation,
               channelCol = "Channel",
               batchCol = "Raw.file",
               suffix = paste0("_TMT", 1:6),
               removeEmptyCols = TRUE)

# check data
head(colData(scp))

#plot rawfile
plot(scp)

# workflow as per 
# https://uclouvain-cbio.github.io/SCP.replication/articles/SCoPE2.html

peptides <- scp[["Peptides"]]
proteins <- scp[["Proteins"]]


#filter out failed runs, as per PSM content
nPSMs <- dims(scp)[1, ]


library(ggplot2)
# plot the number of peptide sequence matches
ggplot(data.frame(nPSMs)) +
  aes(x = nPSMs) +
  geom_histogram() +
  geom_vline(xintercept = 500)

# would drop all assays with a psm lower than 500
scp <- scp[, , nPSMs > 500]

scp <- pep2qvalue(scp,
                  i = names(scp),
                  PEP = "PEP",
                  rowDataName = "qvalue_psm")

scp <- pep2qvalue(scp,
                  i = names(scp),
                  groupBy = "Protein.names",
                  PEP = "PEP",
                  rowDataName = "qvalue_protein")

library(tidyr)

rbindRowData(scp, i = names(scp)) %>%
  data.frame %>%
  pivot_longer(cols = c("PEP", "qvalue_psm", "qvalue_protein"),
               names_to = "measure") %>%
  
  ggplot(aes(x = value)) +
  geom_histogram() +
  geom_vline(xintercept = 0.1) +
  scale_x_log10() +
  facet_grid(rows = vars(measure))

# filter by 1% PSM and protein FDR

scp <- filterFeatures(scp,
                      ~ qvalue_psm < 0.01 & qvalue_protein < 0.01)


#filter out contaminants (all proteins matching to contaminant (starting with CON) or decoy are deleted)
scp <- filterFeatures(scp,
                      ~ !grepl("REV|CON", protein))


#filter noisy spectra (PIF = parental ion fraction, co-isolated peptides)

scp <- filterFeatures(scp,
                      ~ !is.na(PIF) & PIF > 0.8)

# compute SCR (sample to carrier ratio)
scp <- computeSCR(scp,
                  i = names(scp),
                  colvar = "SampleType",
                  carrierPattern = "Carrier",
                  samplePattern = 4:16,
                  rowDataName = "MeanSCR")

# plot SCR rates as histogram

rbindRowData(scp, i = names(scp)) %>%
  data.frame %>%
  ggplot(aes(x = MeanSCR)) +
  geom_histogram() +
  geom_vline(xintercept = c(1/200, 0.1),
             lty = 2:1) +
  scale_x_log10()

scp <- filterFeatures(scp,
                      ~ !is.na(MeanSCR) &
                        !is.infinite(MeanSCR) &
                        MeanSCR < 0.1)

# normalize to reference

scp <- divideByReference(scp,
                         i = names(scp),
                         colvar = "SampleType",
                         samplePattern = ".",
                         refPattern = "Reference")


# aggregate PSM 
remove.duplicates <- function(x)
  apply(x, 2, function(xx) xx[which(!is.na(xx))[1]] )


# give name to the aggregated peptides
peptideAssays <- paste0("peptides_", names(scp))

# aggregate PSMs in different batches to peptides
scp <- aggregateFeaturesOverAssays(scp,
                                   i = names(scp),
                                   fcol = "peptides",
                                   name = peptideAssays,
                                   fun = remove.duplicates)

# join all sets into one assay
# razor peptides (=found in different groups) will be assigned to the group with the most peptides
rbindRowData(scp, i = names(scp)[1:173]) %>%
  data.frame %>%
  group_by(peptide) %>%
  ## The majority vote happens here
  mutate(protein = names(sort(table(protein),
                              decreasing = TRUE))[1]) %>%
  select(peptide, protein) %>%
  filter(!duplicated(peptide, protein)) ->
  ppMap
consensus <- lapply(peptideAssays, function(i) {
  ind <- match(rowData(scp[[i]])$peptide, ppMap$peptide)
  DataFrame(protein = ppMap$protein[ind])
})
names(consensus) <- peptideAssays
rowData(scp) <- consensus


# clean missing data
scp <- infIsNA(scp, i = peptideAssays)
scp <- zeroIsNA(scp, i = peptideAssays)

#join assays
scp <- joinAssays(scp,
                  i = peptideAssays,
                  name = "peptides")

#filter single-cells based on median CV
scp <- medianCVperCell(scp,
                       i = peptideAssays,
                       groupBy = "protein",
                       nobs = 6,
                       na.rm = TRUE,
                       colDataName = "MedianCV",
                       norm = "SCoPE2")







