library(limma)
library(dplyr)

exp_matrix <- read.delim("/home/lukas/Downloads/Raw_data.txt")

design_matrix <- read.delim("/home/lukas/Downloads/Design.txt")

# log transform expression values
exp_matrix[,-1] <- log10(exp_matrix[,-1]) 

# #normalize colwise
exp_matrix_norm <- sweep(exp_matrix[,-1], 2, colSums(exp_matrix[,-1]), FUN="-")
# #normalize rowwise
exp_matrix_norm <- sweep(exp_matrix[,-1], 1, rowSums(exp_matrix[,-1]), FUN="-")
 
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)

corr_matrix <- cor(exp_matrix_norm[,-1])
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
batch <- design_matrix$Batch

library(sva)

combat_edata1 = ComBat(dat=exp_matrix_norm, batch=batch)

# or use limma
y2 <- removeBatchEffect(exp_matrix_norm, batch)

# compare corrected batch effects
boxplot(as.data.frame(exp_matrix_norm),main="Original", ylim=c(0.01,0.05))
boxplot(as.data.frame(y2),main="Batch corrected with Limma", ylim=c(0.01,0.05))
boxplot(as.data.frame(combat_edata1),main="Batch corrected with Combat", ylim=c(0.01,0.05))

# Perform differential expression analysis between each of the two groups

# drop pool columns
exp_matrix_de <- data.frame(combat_edata1) %>% dplyr::select(-starts_with('Control'))
exp_matrix_de <- as.matrix(exp_matrix_de)
rownames(exp_matrix_de) <- (exp_matrix$ID)

# create design matrix consisting of group variables
design_samples <- subset(design_matrix, Group != "Pool")
design <- model.matrix(~0+factor(design_samples$Group))
colnames(design) <- c("EC", "HC", "VP")

# Fit the expression matrix to a linear model
fit <- lmFit(exp_matrix_de, design)

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
exp_matrix_de_paired <- data.frame(exp_matrix_de) %>% dplyr::select(-starts_with('HC'))
design_samples_paired <- design_samples %>%  dplyr::filter(Group!='HC')

# define pairs --> assuming integers indicate pairs
pairs <- factor(c(seq(1,9),seq(1,9)))

# define group
group <- factor(design_samples_paired$Group)

# create design matrix
paired_design <- model.matrix(~pairs + group)

# apply linear model
fit_2 <- lmFit(exp_matrix_de_paired, paired_design)

fit_2 <- eBayes(fit_2)
# shows fold expression and p-values
topTable(fit_2, coef="groupVP")
volcanoplot(fit_2, main="EC versus VP")



fetched_factor <- design_samples_paired[c("Batch", "Gender")]

factors <- lapply(fetched_factor, factor)

factor_character <- paste0("factors['", names(fetched_factor), "']", collapse = "+")

paired_design <- model.matrix(~ factors$Batch + factors$Gender + group)




