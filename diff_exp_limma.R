library(limma)
library(dplyr)

evidence_data <- read.delim("/home/lukas/Downloads/Raw_data.txt")

meta_data_0 <- read.delim("/home/lukas/Downloads/Design.txt")

rownames(evidence_data) <- evidence_data$ID
evidence_data <- evidence_data[ , !(names(evidence_data) %in% c("ID"))]

selectedSampleType_to_exclude <- c("Pool")
evidence_data <- evidence_data[, !(meta_data_0$Group %in%  selectedSampleType_to_exclude)]      
meta_data_0 <- meta_data_0[!(meta_data_0$Group %in%  selectedSampleType_to_exclude), ]


barplot(log10(dim(evidence_data)[1]), main = "log Count of rows")

evidence_data_medians <- data.frame(median_intensity = colMedians(as.matrix(evidence_data)))

evidence_data_medians %>%
  ggplot() +
  aes(x = median_intensity, 
      y = meta_data_0$Group,
      fill = meta_data_0$Group) +
  geom_boxplot() +
  scale_x_log10() +
  labs(fill=as.character(meta_data_0$Group), y=as.character(meta_data_0$Group))





# log transform expression values
evidence_data <- log10(evidence_data) 

protein_matrix <- as.matrix(evidence_data)
protein_matrix <- sweep(protein_matrix, 2, colMedians(protein_matrix), FUN="/")
# #normalize rowwise
protein_matrix <- sweep(protein_matrix, 1, rowMeans(protein_matrix), FUN="/")
evidence_data <- data.frame(protein_matrix)


evidence_data <- replace(evidence_data, is.na(evidence_data), median(unlist(evidence_data), na.rm = TRUE))

library('corrr')
library("ggcorrplot")
library("FactoMineR")
library("factoextra")

corr_matrix <- cor(evidence_data)
ggcorrplot(corr_matrix)

# in this case a correlation within the batches can be clearly observed 

# apply pca on covariance matrix
data.pca <- princomp(corr_matrix)
summary(data.pca)
# scree plot for explained variance by principle component
# in this case the first compononents explain the majority of variance 
fviz_eig(data.pca)

# in the biplot the batch effect can be further evaluated.
# samples within batches show the same vectors
fviz_pca_var(data.pca, col.var = "black")



# get results of pca
fviz_pca_ind(data.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

# batch effects obtained --> use combat to correct
batch <- meta_data_0$Batch

library("scater")
dimred_pca <- calculatePCA(evidence_data,                                          ncomponents = 5,
                           ntop = Inf,
                           scale = TRUE)

dimred_pca <- data.frame(dimred_pca)

ggplot(dimred_pca, aes(x=PC1, y=PC2)) +
  geom_jitter()



library(sva)

combat_edata1 = ComBat(dat=evidence_data_norm, batch=batch)

# or use limma
y2 <- removeBatchEffect(evidence_data_norm, batch)

# compare corrected batch effects
boxplot(as.data.frame(evidence_data_norm),main="Original", ylim=c(0.01,0.05))
boxplot(as.data.frame(y2),main="Batch corrected with Limma", ylim=c(0.01,0.05))
boxplot(as.data.frame(combat_edata1),main="Batch corrected with Combat", ylim=c(0.01,0.05))

# Perform differential expression analysis between each of the two groups

# drop pool columns
evidence_data_de <- data.frame(combat_edata1) %>% dplyr::select(-starts_with('Control'))
evidence_data_de <- as.matrix(evidence_data_de)
rownames(evidence_data_de) <- (evidence_data$ID)

# create design matrix consisting of group variables
design_samples <- subset(meta_data_0, Group != "Pool")
design <- model.matrix(~0+factor(design_samples$Group))
colnames(design) <- c("EC", "HC", "VP")

# Fit the expression matrix to a linear model
fit <- lmFit(evidence_data_de, design)

# EC vs HC
cont_matrix_1 <- makeContrasts(ECvsHC = EC-HC,levels=design)
# Compute contrast
fit_contrast_1 <- contrasts.fit(fit, cont_matrix_1)
# Bayes statistics of differential expression
fit_contrast_1 <- eBayes(fit_contrast_1)
# Generate a vocalno plot to visualize differential expression
volcanoplot(fit_contrast_1, main= "EC versus HC")
# Generate a list of top 100 differentially expressed genes
top_genes_1 <- topTable(fit_contrast_1, number = 100, adjust = "BH")
# Summary of results (number of differentially expressed genes)
result_1 <- decideTests(fit_contrast_1)
summary(result_1)

# EC vs VP
cont_matrix_2 <- makeContrasts(ECvsVP = EC-VP,levels=design)
# Compute contrast
fit_contrast_2 <- contrasts.fit(fit, cont_matrix_2)
# Bayes statistics of differential expression
fit_contrast_2 <- eBayes(fit_contrast_2)
# Generate a vocalno plot to visualize differential expression
volcanoplot(fit_contrast_2, main= "EC versus VP")
# Generate a list of top 100 differentially expressed genes
top_genes_2 <- topTable(fit_contrast_2, number = 100, adjust = "BH")
# Summary of results (number of differentially expressed genes)
result_2 <- decideTests(fit_contrast_2)
summary(result_2)

# HC vs VP
cont_matrix_3 <- makeContrasts(HCvsVP = HC-VP,levels=design)
# Compute contrast
fit_contrast_3 <- contrasts.fit(fit, cont_matrix_3)
# Bayes statistics of differential expression
fit_contrast_3 <- eBayes(fit_contrast_3)
# Generate a vocalno plot to visualize differential expression
volcanoplot(fit_contrast_3, main="HC versus VP")
# Generate a list of top 100 differentially expressed genes
top_genes_3 <- topTable(fit_contrast_3, number = 100, adjust = "BH")
# Summary of results (number of differentially expressed genes)
result_3 <- decideTests(fit_contrast_3)
summary(result_3)


# perform paired group-wise analysis between EC and VP
# exclude HC
evidence_data_de_paired <- data.frame(evidence_data_de) %>% dplyr::select(-starts_with('HC'))
design_samples_paired <- design_samples %>%  dplyr::filter(Group!='HC')

# define pairs --> assuming integers indicate pairs
pairs <- factor(c(seq(1,9),seq(1,9)))

# define group
group <- factor(design_samples_paired$Group)

# create design matrix
paired_design <- model.matrix(~pairs + group)

# apply linear model
fit_2 <- lmFit(evidence_data_de_paired, paired_design)

fit_2 <- eBayes(fit_2)
# shows fold expression and p-values
topTable(fit_2, coef="groupVP")
volcanoplot(fit_2, main="EC versus VP")



fetched_factor <- design_samples_paired[c("Batch", "Gender")]

factors <- lapply(fetched_factor, factor)

factor_character <- paste0("factors['", names(fetched_factor), "']", collapse = "+")

paired_design <- model.matrix(~factors$Gender + group)




