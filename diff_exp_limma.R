library(limma)
library(dplyr)
library(ggplot2)
library(matrixStats)
# load files
evidence_data <- read.delim("/home/lukas/Desktop/MS-Data/example/expression_matrix/Raw_data.txt")

meta_data_0 <- read.delim("/home/lukas/Desktop/MS-Data/example/expression_matrix/Design.txt")


# handle id column
rownames(evidence_data) <- evidence_data$ID
evidence_data <- evidence_data[ , !(names(evidence_data) %in% c("ID"))]




# handle exclusion dependencies
selectedSampleType_to_exclude <- c("Pool")
      
meta_data_0 <- meta_data_0[!(meta_data_0$Group %in%  selectedSampleType_to_exclude), ]


evidence_data <- evidence_data[ , which(meta_data_0$ID == colnames(evidence_data))]
# plot basic counts
#barplot(log10(dim(evidence_data)[1]), main = "log Count of rows")

# # plot median intensities 
# evidence_data_medians <- data.frame(median_intensity = colMedians(as.matrix(evidence_data)))
# 
# evidence_data_medians %>%
#   ggplot() +
#   aes(x = median_intensity, 
#       y = meta_data_0$Group,
#       fill = meta_data_0$Group) +
#   geom_boxplot() +
#   scale_x_log10() +
#   labs(fill=as.character(meta_data_0$Group), y=as.character(meta_data_0$Group))

# plot covariance
#heatmap(cor(t(evidence_data)))

# feature wise output across channels
# selectedProtein <- "P26640"
# to_plot <- data.frame(t(evidence_data[selectedProtein, ]))
# to_plot$names <- rownames(to_plot)
# 
# library(reshape2)
# 
# d = melt(to_plot, id.vars = "names")
# 
# ggplot(data = d,
#        mapping = aes(x = names, y = value, fill=meta_data_0$Batch)) + 
#   geom_col(position = position_dodge()) +
#   labs(y=paste("intensity of", selectedProtein))

# log transform expression values
evidence_data <- log10(evidence_data) 

protein_matrix <- as.matrix(evidence_data)
protein_matrix <- sweep(protein_matrix, 2, colMedians(protein_matrix), FUN="-")
# #normalize rowwise
protein_matrix <- sweep(protein_matrix, 1, rowMeans(protein_matrix), FUN="-")
evidence_data <- data.frame(protein_matrix)


# evidence_data <- replace(evidence_data, is.na(evidence_data), median(unlist(evidence_data), na.rm = TRUE))


# # test MA plot and constand
# 
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
# exp_matrix_0 <- as.matrix(evidence_data)
# colnames(exp_matrix_0) <- meta_data_0$Group
# 
# user_choice <- "EC-HC"
# # split user choice of comp back to sample types
# user_choice_vector <- strsplit(user_choice, split = "-")
# # and assign them to a variable
# choice_A <- user_choice_vector[[1]][1]
# choice_B <- user_choice_vector[[1]][2]
# 
# #find the row indeces of corresponding to the individual sample types
# st_indeces <- split(seq_along(meta_data_0$Group), meta_data_0$Group)
# index_A <- st_indeces[choice_A]
# index_B <- st_indeces[choice_B]
# 
# MAplot(exp_matrix_0[,index_A[[1]]], exp_matrix_0[,index_B[[1]]])








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


# library("scater")
# dimred_pca <- calculatePCA(evidence_data,                                          ncomponents = 5,
#                            ntop = Inf,
#                            scale = TRUE)
# 
# dimred_pca <- data.frame(dimred_pca)
# 
# # sample continuos variable to test 
# 
# meta_data_0$body_weight <- rnorm(nrow(meta_data_0), mean=65)
# 
# ggplot(dimred_pca, aes(x=PC1, 
#                        y=PC2, 
#                        color=meta_data_0$Batch,
#                        shape = meta_data_0$Gender,
#                        size = meta_data_0$body_weight)) +
#   geom_point(alpha=3/4)
# 
# 
# dimred_umap <- calculateUMAP(evidence_data,
#                              ncomponents = 5,
#                              ntop = Inf,
#                              scale = TRUE)
# 
# dimred_umap <- data.frame(dimred_umap)
# 
# ggplot(dimred_umap, aes(x=UMAP1, 
#                         y=UMAP2, 
#                         color=meta_data_0$Batch,
#                         shape = meta_data_0$Gender,
#                         size = meta_data_0$body_weight)) +
#   geom_point(alpha=3/4)
# 
# plot_ly(x=dimred_umap$UMAP1,
#         y=dimred_umap$UMAP2, 
#         z=dimred_umap$UMAP3, 
#         type="scatter3d", 
#         mode="markers",
#         color=meta_data_0$Batch) 


library(sva)

# batch effects obtained --> use combat to correct
batch <- meta_data_0$Batch


combat_edata1 = ComBat(dat=evidence_data, batch=batch)

# or use limma
y2 <- removeBatchEffect(evidence_data, batch)

# compare corrected batch effects
boxplot(as.data.frame(evidence_data),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected with Limma")
boxplot(as.data.frame(combat_edata1),main="Batch corrected with Combat")

# Perform differential expression analysis between each of the two groups
evidence_data_de <- combat_edata1


# create design matrix consisting of group variables
design <- model.matrix(~0+factor(meta_data_0$Group))
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


## Gene set enrichment analysis (GSEA)
# # install packages
# BiocManager::install("fgsea")
# BiocManager::install("clusterProfiler")

library("fgsea")
library("clusterProfiler")
library("edgeR")


# top table of significant genes is the starting point of the analysis
tt <- top_genes_2

mask <- tt$adj.P.Val < 0.05 &
  abs(tt$logFC) > 0.0005
deGenes <- rownames(tt)[mask]
head(deGenes)

geneUniverse <- rownames(tt)
length(geneUniverse)

library(clusterProfiler)

ans.kegg <- enrichKEGG(
  deGenes,
  organism = "hsa",
  keyType = "uniprot",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

tab.kegg <- as.data.frame(ans.kegg)
tab.kegg<- subset(tab.kegg, Count>5)
tab.kegg[1:5, 1:6]

library(enrichplot)
p1 <- graphics::barplot(ans.kegg, showCategory=10)
p1

p2 <- dotplot(ans.kegg, showCategory=20) + ggtitle("KEGG")
p2

cowplot::plot_grid(p1, p3, p5, ncol=2, labels=LETTERS[1:3])