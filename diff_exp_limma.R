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
evidence_data <- log2(evidence_data) 
evidence_data_transf <- evidence_data


protein_matrix <- as.matrix(evidence_data)
protein_matrix <- sweep(protein_matrix, 2, colMedians(protein_matrix), FUN="-")
# #normalize rowwise
protein_matrix <- sweep(protein_matrix, 1, rowMeans(protein_matrix), FUN="-")
evidence_data <- data.frame(protein_matrix)



# # compare mutual info between batches 
# batch_select1 <- "Batch_1"
# batch_select2 <- "Batch_2"
# batch_select3 <- "Batch_3"
# 
# evidence_data_btch1 <- evidence_data[ , which(meta_data_0$Batch == batch_select1)]
# evidence_data_btch1 <- data.frame(rowMeans(evidence_data_btch1))
# 
# evidence_data_btch2 <- evidence_data[ , which(meta_data_0$Batch == batch_select2)]
# evidence_data_btch2 <- data.frame(rowMeans(evidence_data_btch2))
# 
# evidence_data_btch3 <- evidence_data[ , which(meta_data_0$Batch == batch_select3)]
# evidence_data_btch3 <- data.frame(rowMeans(evidence_data_btch3))
# 
# 
# btch_comp_df <- cbind(evidence_data_btch1, evidence_data_btch2, evidence_data_btch3)
# 
# 
# data <- as.matrix(btch_comp_df)
# splits = c(1,1)
# 
# entropy <- lnn_mi(data, splits)
# 
# 

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








# library('corrr')
# library("ggcorrplot")
# library("FactoMineR")
# library("factoextra")
# 
# corr_matrix <- cor(evidence_data)
# ggcorrplot(corr_matrix)
# 
# # in this case a correlation within the batches can be clearly observed 
# 
# # apply pca on covariance matrix
# data.pca <- princomp(corr_matrix)
# summary(data.pca)
# # scree plot for explained variance by principle component
# # in this case the first compononents explain the majority of variance 
# fviz_eig(data.pca)
# 
# # in the biplot the batch effect can be further evaluated.
# # samples within batches show the same vectors
# fviz_pca_var(data.pca, col.var = "black")
# 
# 
# 
# # get results of pca
# fviz_pca_ind(data.pca, col.ind = "cos2", 
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE # Avoid text overlapping (slow if many points)
# )


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

evidence_data <- data.frame(evidence_data)
# batch effects obtained --> use combat to correct
batch <- meta_data_0$Batch
group <- meta_data_0$Group

model <- model.matrix(~group, data = evidence_data)

Combat_batchC <- function(i_exp_matrix, i_batch, i_model) {
  i_exp_matrix <- tryCatch(
    {
      i_exp_matrix <- ComBat(dat = i_exp_matrix,
                             batch = i_batch,
                             mod = i_model)
      return(i_exp_matrix)
    },
    error = function(cond) {
      i_exp_matrix <- ComBat(dat = i_exp_matrix,
                             batch = i_batch)
      print(paste("confounder detected! just batch corrected. reconsider your experimental design", cond))
      return(i_exp_matrix)
    }
  )
}

qr(model)$rank < ncol(model)



combat_edata = Combat_batchC(evidence_data, batch, model)

# # or use limma
y2 <- removeBatchEffect(evidence_data, batch)
# 
# # compare corrected batch effects
# boxplot(as.data.frame(evidence_data),main="Original")
# boxplot(as.data.frame(y2),main="Batch corrected with Limma")
# boxplot(as.data.frame(combat_edata1),main="Batch corrected with Combat")

# Perform differential expression analysis between each of the two groups
evidence_data <- combat_edata


# create design matrix consisting of group variables
design <- model.matrix(~0+factor(meta_data_0$Group) + factor(meta_data_0$Batch))
colnames(design) <- c("EC", "HC", "Pool", "VP", "Batch_2", "Batch_3")


# Fit the expression matrix to a linear model
fit <- lmFit(evidence_data, design)
bayes_fit <- eBayes(fit)


# EC vs HC
cont_matrix_1 <- makeContrasts(contrasts ="EC-HC",levels=design)
# Compute contrast
fit_contrast_1 <- contrasts.fit(fit, cont_matrix_1)
# Bayes statistics of differential expression
fit_contrast_1 <- eBayes(fit_contrast_1)
# Generate a vocalno plot to visualize differential expression
volcanoplot(fit_contrast_1, main= "EC versus HC")
# Generate a list of top 100 differentially expressed genes
top_genes_1 <- topTable(fit_contrast_1, number = Inf, adjust = "BH")
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
top_genes_2 <- topTable(fit_contrast_2, number = Inf, adjust = "BH")
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
top_genes_3 <- topTable(fit_contrast_3, number = Inf, adjust = "BH")
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
tt <- top_genes_1

background_genes <- data.frame(rownames(tt))

library(tidyverse)
write_delim(
  background_genes,
  "/home/lukas/Downloads/background_genes.txt",
  delim = " ",
)

background <- read.delim("/home/lukas/Downloads/background_genes.txt") 
background <- as.vector(background[1])
background <- background[[1]]

mask <- tt$adj.P.Val < 1 &
  abs(tt$logFC) > 0.0005
deGenes <- rownames(tt)[mask]
head(deGenes)

ans.kegg <- enrichKEGG(
  deGenes,
  organism = "hsa",
  keyType = "uniprot",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE,
  universe = background
)

tab.kegg <- as.data.frame(ans.kegg)
#tab.kegg<- subset(tab.kegg, Count>5)
# 
# tab.kegg$GeneRatio[1] <- "10/53"
# tab.kegg$Count[1] <- 10

tab.kegg$denominator <- as.numeric(gsub("^\\d+/(\\d+)$", "\\1", tab.kegg$GeneRatio))
tab.kegg$decimal_gene_ratio <- tab.kegg$Count / tab.kegg$denominator

library(enrichplot)
p1 <- graphics::barplot(ans.kegg, showCategory=10)
p1

p2 <- dotplot(ans.kegg, showCategory=20) + ggtitle("KEGG")
p2

cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])

library(plotly)
p3 <- plot_ly(data=tab.kegg, 
              x=~decimal_gene_ratio, 
              y=~Description, 
              type = "scatter", 
              color= ~p.adjust, 
              size=~decimal_gene_ratio, 
              text=~GeneRatio,
              hovertemplate= paste('%{y}', '<br>Gene ratio: %{text}<br><extra></extra>')) %>%
  layout(xaxis=list(
    title="Gene Ratio"
    ))
  
p3


library("org.Hs.eg.db")
library("AnnotationDbi")

myPval <- data.frame(p.val=tt$P.Value, t.val=tt$t)
myPval$UNIPROT <- rownames(tt)

entrez_ids <- select(org.Hs.eg.db, myPval$UNIPROT, "ENTREZID", "UNIPROT")

myPval <- merge(myPval, entrez_ids, by="UNIPROT")

myPval2 <- data.frame(myPval$t.val)
myPval2$ENTREZID <- myPval$ENTREZID

GO_enrichment <- enrichGO(
  myPval2$ENTREZID,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)

plotGOgraph(GO_enrichment)









library("piano")

# parse data
myPval <- data.frame(p.val=tt$P.Value, t.val=tt$t, logFC=tt$logFC)
myPval$UNIPROT <- rownames(tt)

library("org.Hs.eg.db")
library("AnnotationDbi")

# fetch Entrez IDs
entrez_ids <- select(org.Hs.eg.db, myPval$UNIPROT, "ENTREZID", "UNIPROT")
myPval <- merge(myPval, entrez_ids, by="UNIPROT")

# check nNA
nrow(myPval)
nrow(myPval[which(is.na(myPval$ENTREZID)) , ])

# check duplicates 
nrow(myPval)
length(unique(myPval$ENTREZID))

myPval <- myPval[which(!is.na(myPval$ENTREZID)) , ]
myPval <- myPval[!duplicated(myPval$ENTREZID) , ]


rownames(myPval) <- myPval$ENTREZID

# gene set collection
myGSC <- loadGSC("/home/lukas/Downloads/c7.all.v2023.1.Hs.entrez.gmt")


logFCs <- data.frame(myPval$logFC)
rownames(logFCs) <- rownames(myPval)

pVals <- data.frame(myPval$p.val)
rownames(pVals) <- rownames(myPval)

library("snowfall")
library("parallel")

cores <- detectCores()

gsaRes <- runGSA(geneLevelStats = pVals,
                 directions = logFCs,
                 gsc=myGSC, 
                 ncpus = cores)


gsa_results <- GSAsummaryTable(gsaRes = gsaRes)


par(mar = c(1, 1, 1, 1))
networkPlot(gsaRes,class="non")


nw <- networkPlot2(gsaRes, class="non", significance = 0.5, shiny = T)

library("visNetwork")

visNetwork(nodes = nw$x$nodes, edges = nw$x$edges)






cut <- nrow(tt[tt$adj.P.Value < 1 & abs(tt$logFC) > 0.0001])

GSAheatmap(gsaRes, cutoff=cut)




