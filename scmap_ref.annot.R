set.seed(9742)
#devtools::install_github("hemberg-lab/scmap")
library(SingleCellExperiment)
library(scmap)
library(tidyverse)
library(scater)
library(Seurat)
library(patchwork)

## Create an `SingleCellExperiment` object containing reference data

data.f3 = Read10X(data.dir = "Resources/filtered_gene_bc_matrices_mex_Fig3/mm10/")
data.f3 <- CreateSeuratObject(counts = data.f3, project = "figure3", 
                              min.cells = 10, min.features = 200) ## dim 14578 18657
data.f3[["percent.mt"]] <- PercentageFeatureSet(data.f3, pattern = "^mt-")
data.f3$log10GenesPerUMI <- log10(data.f3$nFeature_RNA) / log10(data.f3$nCount_RNA) # 998.3 MB
data.f3 <- subset(data.f3, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8)
data.f3 <- NormalizeData(data.f3, normalization.method = "LogNormalize", scale.factor = 10000) # 989.9 MB

meta = read.csv("Resources/References/annot_Fig_3.csv") %>% dplyr::select(cell, cluster)
str = colnames(data.f3)
str = str[str %in% meta$cell]
data.f3 = data.f3[,str]
meta = meta[meta$cell %in% str, ]
rownames(meta) <- meta$cell
meta$cell <- NULL

sce = SingleCellExperiment(assays = list(logcounts = as.matrix(assay(as.SingleCellExperiment(data.f3), "logcounts"))), 
                           colData = meta)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
sce

############## 0901 ver. #############

data.f3.mtx = as.matrix(GetAssayData(data.f3))
sce = SingleCellExperiment(data.f3, colData = meta)
meta = read.csv("Resources/References/annot_Fig_3.csv") %>% dplyr::select(cell, cluster)
str = colnames(data.f3.mtx)
str = str[str %in% meta$cell] 
data.f3.mtx = data.f3.mtx[,str]

meta = meta[meta$cell %in% str, ]
rownames(meta) <- meta$cell
meta$cell <- NULL

sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(data.f3.mtx)), colData = meta)
logcounts(sce) <- log2(normcounts(sce) + 1)
# use gene names as feature symbols
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rownames(sce)), ]
sce

#######################################

sce <- selectFeatures(sce, suppress_plot = FALSE)
sce <- indexCluster(sce, cluster_col = "cluster")
scmap_cluster_reference <- metadata(sce)$scmap_cluster_index
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))

saveRDS(sce, "data.f3.reference_scmap.0902.RData")
sce = readRDS("data.f3.reference_scmap.0902.RData")

harmony = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")
har.sce <- as.SingleCellExperiment(harmony)
# use gene names as feature symbols
rowData(har.sce)$feature_symbol <- rownames(har.sce)


## Select the most informative features (genes) from input dataset
har.sce <- selectFeatures(har.sce, suppress_plot = FALSE) # Features highlighted with the red colour will be used in the futher analysis (projection).
table(rowData(har.sce)$scmap_features)
har.sce <- indexCluster(har.sce, cluster_col = "seurat_clusters")
head(metadata(har.sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(har.sce)$scmap_cluster_index))

scmapCluster_results <- scmapCluster(
  projection = har.sce, 
  index_list = list(immune1 = scmap_cluster_reference), 
  threshold=0.1
  )
# plot the results of our annotation
par(mar=c(13, 4, 1, 0))
barplot(table(scmapCluster_results$combined_labs), las=2)

# Store this annotation information within the query object
harmony@meta.data$scmap_cluster <- scmapCluster_results$combined_labs

# Make a UMAP of the cells, labeled with the cell-type annotations from scmapCluster
DimPlot(harmony, reduction = "umap", group.by = "scmap_cluster",label = T, label.size = 5, repel = 0.5)
scater::plotReducedDim(harmony, dimred="UMAP", colour_by="scmap_cluster")






## Create an `SingleCellExperiment` object containing reference data

data.f3 = Read10X(data.dir = "Resources/filtered_gene_bc_matrices_mex_")
data.f3 <- CreateSeuratObject(counts = data.f3, project = "figure3", 
                              min.cells = 10, min.features = 200) ## dim 14578 18657
data.f3[["percent.mt"]] <- PercentageFeatureSet(data.f3, pattern = "^mt-")
data.f3$log10GenesPerUMI <- log10(data.f3$nFeature_RNA) / log10(data.f3$nCount_RNA) # 998.3 MB
data.f3 <- subset(data.f3, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8)
data.f3 <- NormalizeData(data.f3, normalization.method = "LogNormalize", scale.factor = 10000) # 989.9 MB

meta = read.csv("Resources/References/annot_Fig_3.csv") %>% dplyr::select(cell, cluster)
rownames(meta) = meta$cell
meta$cell <- NULL
str = colnames(data.f3)
str = str[str %in% rownames(meta)]
data.f3 = data.f3[,str]
meta = meta[rownames(meta) %in% str, ]

sce = SingleCellExperiment(assays = list(logcounts = as.matrix(assay(as.SingleCellExperiment(data.f3), "logcounts"))), 
                           colData = meta)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
sce

############## 0901 ver. #############

data.f3.mtx = as.matrix(GetAssayData(data.f3))
sce = SingleCellExperiment(data.f3, colData = meta)
meta = read.csv("Resources/References/annot_Fig_3.csv") %>% dplyr::select(cell, cluster)
str = colnames(data.f3.mtx)
str = str[str %in% meta$cell] 
data.f3.mtx = data.f3.mtx[,str]

meta = meta[meta$cell %in% str, ]
rownames(meta) <- meta$cell
meta$cell <- NULL

sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(data.f3.mtx)), colData = meta)
logcounts(sce) <- log2(normcounts(sce) + 1)
# use gene names as feature symbols
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rownames(sce)), ]
sce

#######################################

sce <- selectFeatures(sce, suppress_plot = FALSE)
sce <- indexCluster(sce, cluster_col = "cluster")
scmap_cluster_reference <- metadata(sce)$scmap_cluster_index
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))

saveRDS(sce, "data.f3.reference_scmap.0902.RData")
sce = readRDS("data.f3.reference_scmap.0902.RData")

harmony = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")
har.sce <- as.SingleCellExperiment(harmony)
# use gene names as feature symbols
rowData(har.sce)$feature_symbol <- rownames(har.sce)


## Select the most informative features (genes) from input dataset
har.sce <- selectFeatures(har.sce, suppress_plot = FALSE) # Features highlighted with the red colour will be used in the futher analysis (projection).
table(rowData(har.sce)$scmap_features)
har.sce <- indexCluster(har.sce, cluster_col = "seurat_clusters")
head(metadata(har.sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(har.sce)$scmap_cluster_index))

scmapCluster_results <- scmapCluster(
  projection = har.sce, 
  index_list = list(immune1 = scmap_cluster_reference), 
  threshold=0.1
)
# plot the results of our annotation
par(mar=c(13, 4, 1, 0))
barplot(table(scmapCluster_results$combined_labs), las=2)

# Store this annotation information within the query object
harmony@meta.data$scmap_cluster <- scmapCluster_results$combined_labs

# Make a UMAP of the cells, labeled with the cell-type annotations from scmapCluster
DimPlot(harmony, reduction = "umap", group.by = "scmap_cluster",label = T, label.size = 5, repel = 0.5)
scater::plotReducedDim(harmony, dimred="UMAP", colour_by="scmap_cluster")



