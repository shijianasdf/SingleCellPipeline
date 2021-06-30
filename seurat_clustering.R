# Tag & sample name : help, project name : Chromium
options(stringsAsFactors = F)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
setwd("~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/")
      
## Load CellRanger outputs  & Setup the Seurat Object
helpless.data <- Read10X(data.dir = "./Data/outs_cellranger_exp/033_filtered_feature_bc_matrix/")
helpless <- CreateSeuratObject(counts = helpless.data, project = "chromium033", 
                               min.cells = 10, min.features = 200) # dim(helples) = 14540 7429
helpless[["percent.mt"]] <- PercentageFeatureSet(helpless, pattern = "^mt-")
helpless$log10GenesPerUMI <- log10(helpless$nFeature_RNA) / log10(helpless$nCount_RNA)

coghelp.data <- Read10X(data.dir = "./Data/outs_cellranger_exp/034_filtered_feature_bc_matrix/")
coghelp <- CreateSeuratObject(counts = coghelp.data, project = "chromium034", 
                              min.cells = 10, min.features = 200) # dim(coghelp) = 14696 10237
coghelp[["percent.mt"]] <- PercentageFeatureSet(coghelp, pattern = "^mt-") 
coghelp$log10GenesPerUMI <- log10(coghelp$nFeature_RNA) / log10(coghelp$nCount_RNA)

sephelp.data <- Read10X(data.dir = "./Data/outs_cellranger_exp/035_filtered_feature_bc_matrix/")
sephelp <- CreateSeuratObject(counts = sephelp.data, project = "chromium035", 
                              min.cells = 10, min.features = 200) #dim(sephelp = 14540 7429)
sephelp[["percent.mt"]] <- PercentageFeatureSet(sephelp, pattern = "^mt-")
sephelp$log10GenesPerUMI <- log10(sephelp$nFeature_RNA) / log10(sephelp$nCount_RNA)

naive.data <- Read10X(data.dir = "./Data/outs_cellranger_exp/040_filtered_feature_bc_matrix/")
naive <- CreateSeuratObject(counts = naive.data, project = "chromium040",
                            min.cells = 10, min.features = 200) # dim(naive) = 14851 8563
naive[["percent.mt"]] <- PercentageFeatureSet(naive, pattern = "^mt-")
naive$log10GenesPerUMI <- log10(naive$nFeature_RNA) / log10(naive$nCount_RNA)



## Merge all data to visualize violinplot and compare
alldata <- merge(helpless, c(coghelp, sephelp, naive), add.cell.ids=c("Chromium033","Chromium034","Chromium035","Chromium040"))
p <- VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot("./vlnplot_alldata_0402.pdf", p, base_height = 6, base_width = 12)

## For Slide
VlnPlot(alldata, features = "percent.mt")

## Pre-processing (QC) and Normalizing Data
helpless <- subset(helpless, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8)
helpless <- NormalizeData(helpless, normalization.method = "LogNormalize", scale.factor = 10000)
coghelp <- subset(coghelp, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8)
coghelp <- NormalizeData(coghelp, normalization.method = "LogNormalize", scale.factor = 10000)
sephelp <- subset(sephelp, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8)
sephelp <- NormalizeData(sephelp, normalization.method = "LogNormalize", scale.factor = 10000)
naive <- subset(naive, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
naive <- NormalizeData(naive, normalization.method = "LogNormalize", scale.factor = 10000)



################################# HELPLESS (Chromium033) ##########################################

## Remove Malat1
helpless = helpless[!row.names(helpless) %in% c("Malat1")] # 14540 5676 -> 14539 5676

## Feature Selection
FeatureScatter(helpless, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
helpless <- FindVariableFeatures(helpless, selection.method = "vst", nfeatures = 5000)
top10 <- head(VariableFeatures(helpless), 10) # Identify the 10 most highly variable genes
plot1 <- VariableFeaturePlot(helpless) # plot variable features 
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

## Scalining the data - linear transformation
all.genes <- rownames(helpless)
helpless <- ScaleData(helpless, features = all.genes)

## Linear dimensional reduction
helpless <- RunPCA(helpless, features = VariableFeatures(object = helpless))
###  several useful ways of visualizing both cells and features that define the PCA
VizDimLoadings(helpless, dims = 1:20, reduction = "pca")
DimPlot(helpless, reduction = "pca")
DimHeatmap(helpless, dims = 1:20, cells = 500, balanced = TRUE)

## Determine the dimensionality of dataset
helpless <- JackStraw(helpless, num.replicate = 100)
helpless <- ScoreJackStraw(helpless, dims = 1:20)
JackStrawPlot(helpless, dims = 1:20) ## ambigious... 
ElbowPlot(helpless) # alternative heuristic method of determining dimensionality.

########## PC 10 ##########
## Cluster the cells 
helpless <- FindNeighbors(helpless, dims = 1:10)
helpless <- FindClusters(helpless, resolution = 0.5) ## increased values -> greater number of clusters

## Run non-linear dimensional reduction (UMAP/tSNE)
helpless <- RunUMAP(helpless, dims = 1:10)
p <- DimPlot(helpless, reduction = "umap", label = T, label.size = 5)
saveRDS(helpless, file = "Data/SeuratRDS/chromium033_PC10_res0.5.rds")
save_plot(paste('Figures/plot_umap', 'PC10', 'helpless', 'png', sep='.'), p, base_height = 5, base_width = 6)

########## PC 10 - tSNE ##########
# ## Cluster the cells
# helpless <- FindNeighbors(helpless, dims = 1:10)
# helpless <- FindClusters(helpless, resolution = 0.5) ## increased values -> greater number of clusters
# 
# ## Run non-linear dimensional reduction (UMAP/tSNE)
# helpless <- RunTSNE(helpless, dims = 1:10)
# p <- DimPlot(helpless, reduction = "tsne")
# #saveRDS(helpless, file = "./chromium033_PC10_res0.5.rds")
# save_plot(paste('Figures/plot_tsne', 'PC10', 'helpless', 'png', sep='.'), p, base_height = 5, base_width = 6)


# find markers for every cluster compared to all remaining cells, report only the positive ones
helpless.markers <- FindAllMarkers(helpless, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(helpless.markers, paste('Tables/MarkerGenes_Seurat/table.FindAllMarkers', 'PC10', 'helpless', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)
t = helpless.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
write.table(t, paste('Tables/MarkerGenes_Seurat/table.Top15FindAllMarkers', 'PC10', 'helpless', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)




################################# COGHELP (Chromium034) ##########################################

## Remove Malat1
coghelp = coghelp[!row.names(coghelp) %in% c("Malat1")] 

## Feature Selection

FeatureScatter(coghelp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
coghelp <- FindVariableFeatures(coghelp, selection.method = "vst", nfeatures = 5000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(coghelp), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(coghelp)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

## Scalining the data - linear transformation
all.genes <- rownames(coghelp)
coghelp <- ScaleData(coghelp, features = all.genes)

## Linear dimensional reduction
coghelp <- RunPCA(coghelp, features = VariableFeatures(object = coghelp))
###  several useful ways of visualizing both cells and features that define the PCA
VizDimLoadings(coghelp, dims = 1:20, reduction = "pca")
DimPlot(coghelp, reduction = "pca")
DimHeatmap(coghelp, dims = 1:20, cells = 500, balanced = TRUE)

## Determine the dimensionality of dataset
coghelp <- JackStraw(coghelp, num.replicate = 100)
coghelp <- ScoreJackStraw(coghelp, dims = 1:20)
JackStrawPlot(coghelp, dims = 1:20) ## ambigious... Let's try 20, 15, 10, 9, 8, 7!
ElbowPlot(coghelp) # alternative heuristic method of determining dimensionality.


########## PC 10 ##########
## Cluster the cells 
coghelp <- FindNeighbors(coghelp, dims = 1:10)
coghelp <- FindClusters(coghelp, resolution = 0.5) ## increased values -> greater number of clusters

## Run non-linear dimensional reduction (UMAP/tSNE)
coghelp <- RunUMAP(coghelp, dims = 1:10)
p <- DimPlot(coghelp, reduction = "umap", label = T, label.size = 5)
saveRDS(coghelp, file = "Data/SeuratRDS/chromium034_PC10_res0.5.rds")
save_plot(paste('Figures/plot_umap', 'PC10', 'coghelp', 'png', sep='.'), p, base_height = 5, base_width = 6)

# find markers for every cluster compared to all remaining cells, report only the positive ones
coghelp.markers <- FindAllMarkers(coghelp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(coghelp.markers, paste('Tables/MarkerGenes_Seurat/table.FindAllMarkers', 'PC10', 'coghelp', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)
t2 = coghelp.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
write.table(t2, paste('Tables/MarkerGenes_Seurat/table.Top15FindAllMarkers', 'PC10', 'coghelp', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)




################################# SEPHELP (Chromium035) ##########################################

## Remove Malat1
sephelp = sephelp[!row.names(sephelp) %in% c("Malat1")] 

## Feature Selection

FeatureScatter(sephelp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sephelp <- FindVariableFeatures(sephelp, selection.method = "vst", nfeatures = 5000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sephelp), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sephelp)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

## Scalining the data - linear transformation
all.genes <- rownames(sephelp)
sephelp <- ScaleData(sephelp, features = all.genes)

## Linear dimensional reduction
sephelp <- RunPCA(sephelp, features = VariableFeatures(object = sephelp))
###  several useful ways of visualizing both cells and features that define the PCA
VizDimLoadings(sephelp, dims = 1:20, reduction = "pca")
DimPlot(sephelp, reduction = "pca")
DimHeatmap(sephelp, dims = 1:20, cells = 500, balanced = TRUE)

## Determine the dimensionality of dataset
sephelp <- JackStraw(sephelp, num.replicate = 100)
sephelp <- ScoreJackStraw(sephelp, dims = 1:20)
JackStrawPlot(sephelp, dims = 1:20) ## ambigious... Let's try 20, 15, 10, 9, 8, 7!
ElbowPlot(sephelp) # alternative heuristic method of determining dimensionality.


########## PC 10 ##########
## Cluster the cells 
sephelp <- FindNeighbors(sephelp, dims = 1:10)
sephelp <- FindClusters(sephelp, resolution = 0.5) ## increased values -> greater number of clusters

## Run non-linear dimensional reduction (UMAP/tSNE)
sephelp <- RunUMAP(sephelp, dims = 1:10)
p <- DimPlot(sephelp, reduction = "umap", label = T, label.size = 5)
saveRDS(sephelp, file = "Data/SeuratRDS/chromium035_PC10_res0.5.rds")
save_plot(paste('Figures/plot_umap', 'PC10', 'sephelp', 'png', sep='.'), p, base_height = 5, base_width = 6)

# find markers for every cluster compared to all remaining cells, report only the positive ones
sephelp.markers <- FindAllMarkers(sephelp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(sephelp.markers, paste('Tables/MarkerGenes_Seurat/table.FindAllMarkers', 'PC10', 'sephelp', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)
t3 = sephelp.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
write.table(t3, paste('Tables/MarkerGenes_Seurat/table.Top15FindAllMarkers', 'PC10', 'sephelp', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)



################################# NAIVE (Chromium040) ##########################################

## Remove Malat1
naive = naive[!row.names(naive) %in% c("Malat1")] 

## Feature Selection
FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
naive <- FindVariableFeatures(naive, selection.method = "vst", nfeatures = 5000)
top10 <- head(VariableFeatures(naive), 10) # Identify the 10 most highly variable genes
plot1 <- VariableFeaturePlot(naive) # plot variable features 
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

## Scalining the data - linear transformation
all.genes <- rownames(naive)
naive <- ScaleData(naive, features = all.genes)

## Linear dimensional reduction
naive <- RunPCA(naive, features = VariableFeatures(object = naive))
###  several useful ways of visualizing both cells and features that define the PCA
VizDimLoadings(naive, dims = 1:20, reduction = "pca")
DimPlot(naive, reduction = "pca")
DimHeatmap(naive, dims = 1:20, cells = 500, balanced = TRUE)

## Determine the dimensionality of dataset
naive <- JackStraw(naive, num.replicate = 100)
naive <- ScoreJackStraw(naive, dims = 1:20)
JackStrawPlot(naive, dims = 1:20) ## ambigious... 
ElbowPlot(naive) # alternative heuristic method of determining dimensionality.



########## PC 10 ##########
## Cluster the cells 
naive <- FindNeighbors(naive, dims = 1:10)
naive <- FindClusters(naive, resolution = 0.5) ## increased values -> greater number of clusters

## Run non-linear dimensional reduction (UMAP/tSNE)
naive <- RunUMAP(naive, dims = 1:10)
p <- DimPlot(naive, reduction = "umap", label = T, label.size = 5)
saveRDS(naive, file = "Data/SeuratRDS/chromium040_PC10_res0.5.rds")
save_plot(paste('Figures/plot_umap', 'PC10', 'naive', 'png', sep='.'), p, base_height = 5, base_width = 6)

# ########## PC 10 - tSNE ##########
# ## Cluster the cells 
# naive <- FindNeighbors(naive, dims = 1:10)
# naive <- FindClusters(naive, resolution = 0.5) ## increased values -> greater number of clusters
# 
# ## Run non-linear dimensional reduction (UMAP/tSNE)
# naive <- RunTSNE(naive, dims = 1:10)
# p <- DimPlot(naive, reduction = "tsne")
# #saveRDS(naive, file = "./chromium033_PC10_res0.5.rds")
# save_plot(paste('Figures/plot_tsne', 'PC10', 'naive', 'png', sep='.'), p, base_height = 5, base_width = 6)

# find markers for every cluster compared to all remaining cells, report only the positive ones
naive.markers <- FindAllMarkers(naive, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(naive.markers, paste('Tables/MarkerGenes_Seurat/table.FindAllMarkers', 'PC10', 'naive', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)
t = naive.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
write.table(t, paste('Tables/MarkerGenes_Seurat/table.Top15FindAllMarkers', 'PC10', 'naive', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)

