options(stringsAsFactors = F)
library(DESeq2)
library(tidyverse)
library(harmony)
library(dplyr)
library(ggplot2)
library(harmony)
library(cowplot)
library(Seurat)
setwd("~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/")


## Harmony vignette 

## Load 3 SeuratObjects
helpless.data <- Read10X(data.dir = "./Data/outs_cellranger_exp/033_filtered_feature_bc_matrix/")
helpless <- CreateSeuratObject(counts = helpless.data, project = "chromium033", 
                               min.cells = 10, min.features = 200) # dim(helples) = 14540 7429
helpless[["percent.mt"]] <- PercentageFeatureSet(helpless, pattern = "^mt-")
helpless$log10GenesPerUMI <- log10(helpless$nFeature_RNA) / log10(helpless$nCount_RNA)
helpless <- subset(helpless, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) # 14540 5676
helpless = helpless[!row.names(helpless) %in% c("Malat1")] # 14540 5676 -> 14539 5676

coghelp.data <- Read10X(data.dir = "./Data/outs_cellranger_exp/034_filtered_feature_bc_matrix/")
coghelp <- CreateSeuratObject(counts = coghelp.data, project = "chromium034", 
                              min.cells = 10, min.features = 200) # dim(coghelp) = 14696 10237
coghelp[["percent.mt"]] <- PercentageFeatureSet(coghelp, pattern = "^mt-") 
coghelp$log10GenesPerUMI <- log10(coghelp$nFeature_RNA) / log10(coghelp$nCount_RNA)
coghelp <- subset(coghelp, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) # 14696 8893
coghelp = coghelp[!row.names(coghelp) %in% c("Malat1")] # 14540 5676 -> 14539 5676

sephelp.data <- Read10X(data.dir = "./Data/outs_cellranger_exp/035_filtered_feature_bc_matrix/")
sephelp <- CreateSeuratObject(counts = sephelp.data, project = "chromium035", 
                              min.cells = 10, min.features = 200) # dim(sephelp) = 14540 7429
sephelp[["percent.mt"]] <- PercentageFeatureSet(sephelp, pattern = "^mt-")
sephelp$log10GenesPerUMI <- log10(sephelp$nFeature_RNA) / log10(sephelp$nCount_RNA)
sephelp <- subset(sephelp, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) # 14842 5776
sephelp = sephelp[!row.names(sephelp) %in% c("Malat1")] # 14540 5676 -> 14539 5676

naive.data <- Read10X(data.dir = "./Data/outs_cellranger_exp/040_filtered_feature_bc_matrix/")
naive <- CreateSeuratObject(counts = naive.data, project = "chromium040",
                            min.cells = 10, min.features = 200) # dim(naive) = 14851 8563
naive[["percent.mt"]] <- PercentageFeatureSet(naive, pattern = "^mt-")
naive$log10GenesPerUMI <- log10(naive$nFeature_RNA) / log10(naive$nCount_RNA)
naive <- subset(naive, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) # 14851 6405
naive = naive[!row.names(naive) %in% c("Malat1")] # 14540 5676 -> 14539 5676

######

set.seed(42)  # seed를 매번 샘플링 할때마다 설정해야 같은 random number (난수)가 발생해서 재현할 수 있습니다.
n_cells = 5000

helpless = subset(helpless, cells = sample(Cells(helpless), n_cells)) # sample 이라는 function으로 sampling을 함.
print (dim(helpless)) # 14539 5000
keep_genes = rownames(helpless)[rowSums(helpless) >= 10]  # sampling하고 나서 유전자 발현량이 0이 되는 유전자 찾기
helpless = subset(helpless, features = keep_genes) #feature filtering
print (dim(helpless)) # 14158 5000
helpless = RenameCells(helpless, add.cell.id = "chromium033")

coghelp = subset(coghelp, cells = sample(Cells(coghelp), n_cells)) #cell sampling
print (dim(coghelp)) # 14695 5000
keep_genes = rownames(coghelp)[rowSums(coghelp) >= 10]
coghelp = subset(coghelp, features = keep_genes) #feature filtering
print (dim(coghelp)) # 13704 5000
coghelp = RenameCells(coghelp, add.cell.id = "chromium034")

sephelp = subset(sephelp, cells = sample(Cells(sephelp), n_cells)) #cell sampling
print (dim(sephelp)) # 14841 5000
keep_genes = rownames(sephelp)[rowSums(sephelp) >= 10]
sephelp = subset(sephelp, features = keep_genes) #feature filtering
print (dim(sephelp)) # 14385 5000
sephelp = RenameCells(sephelp, add.cell.id = "chromium035")

naive = subset(naive, cells = sample(Cells(naive), n_cells)) #cell sampling
print (dim(naive)) # 14850 5000
keep_genes = rownames(naive)[rowSums(naive) >= 10]
naive = subset(naive, features = keep_genes) #feature filtering
print (dim(naive)) # 14227 5000
naive = RenameCells(naive, add.cell.id = "chromium040")

## Make a list and put these three objects - maybe need sampling. 
d.sc = list()
d.sc[['helpless']] <- helpless
d.sc[['coghelp']] <- coghelp
d.sc[['sephelp']] <- sephelp
d.sc[['naive']] <- naive

#saveRDS(d.sc, 'Data/Seuratobjectslist_merged_after_sampling_naive_added_0401.rds') ## Removed Malat1
d.sc = readRDS('Data/Seuratobjectslist_merged_after_sampling_naive_added_0401.rds')


############# Without Naive ###############

sc <- merge(d.sc[["helpless"]], y = unlist(c(d.sc[["coghelp"]], d.sc[["sephelp"]]))) # 14740 15000. Malat1 removed.
sc@meta.data$orig.ident <- factor(sc@meta.data$orig.ident, levels = c("chromium033", "chromium034", "chromium035"))
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 5000)
## Run ScaleData code and save the RDS data - remove after RunHarmony.
sc <- ScaleData(sc, features = rownames(sc), vars.to.regress = c("nFeature_RNA", "nCount_RNA")) ## nFeature,nCount가 다름.
sc <- RunPCA(sc, features = VariableFeatures(object = sc)) 

#saveRDS(sc, paste('Data/mergedsc_runPCA_completed_without_naive_0406.rds', sep='')) 
#sc = readRDS('Data/mergedsc_runPCA_completed_without_naive_0406.rds')

## Figure 
p = VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = 'orig.ident', ncol = 3, pt.size = 0.001)
#save_plot(paste('./Figures/vlnplot.QC_mergedsc_without_naive_0406.png',sep='.'), p, base_height = 7.5, base_width = 15)
plot1 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", pt.size = 0.5)
plot2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p = plot1 + plot2
#save_plot(paste('./Figures/scatterplot.QC_mergedsc_without_naive_0406.png',sep='.'), p, base_height = 6, base_width = 10)

## Run Harmony
options(repr.plot.height = 3, repr.plot.width = 5)
sc1 <- RunHarmony(sc, "orig.ident", plot_convergence = T) # 4 iterations

# Compare the Plot - before & after harmony
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = sc, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sc, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
#save_plot(paste('./Figures/dim.vimplot_mergedsc_without_naive_0406.png',sep='.'), last_plot(), base_height = 5, base_width = 12)
p1 <- DimPlot(object = sc1, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sc1, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
save_plot(paste('./Figures/dim.vimplot_Harmony_without_naive_0406.png',sep='.'), last_plot(), base_height = 5, base_width = 12)

## Perplexity : determine clustering strength (higher value -> more distinct clusters)
ndims = 10
sc1 <- RunTSNE(sc1, reduction = "harmony", dims = 1:ndims, min.dist = 0.1, perplexity = 100) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:ndims) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

#saveRDS(sc1, paste('Harmony_cluster_0to16_per100_without_naive_0406.rds', sep=''))
sc1 <- readRDS('Data/Harmony_cluster_0to16_per100_without_naive_0406.rds')

DimPlot(sc1, reduction = "tsne", label = TRUE, pt.size = .1)
#save_plot(paste('./Figures/TSNE.Harmony_without_naive_all_0406.png',sep='.'), last_plot(), base_height = 6, base_width = 7)
DimPlot(sc1, reduction = "tsne", label = TRUE, pt.size = .1, split.by = "orig.ident")
#save_plot(paste('./Figures/TSNE.Harmony_without_naive_split_0406.png',sep='.'), last_plot(), base_height =5, base_width = 15)

sc1 <- BuildClusterTree(sc1, reorder.numeric = T)
Tool(object = sc1, slot = 'BuildClusterTree')
PlotClusterTree(sc1)
save_plot(paste('./Figures/clustertree.Harmony_without_naive_0406.png',sep='.'), PlotClusterTree(sc1), base_height = 7, base_width = 7)

prop.table(table(Idents(sc1)))
table(Idents(sc1), sc1$orig.ident)
d = as.data.frame(prop.table(table(Idents(sc1), sc1$orig.ident)))
#write.table(d, 'Tables/table.frequecy.Harmony.without_naive.0406.txt', sep='\t', quote = F, row.names = F, col.names = T)

# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(sc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(markers, paste('Tables/table.FindAllMarkers_without_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'), 
#            sep='\t', quote = F, row.names = F, col.names = T)
t = markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
#write.table(t, paste('Tables/table.Top15FindAllMarkers_without_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'), 
#            sep='\t', quote = F, row.names = F, col.names = T)





############# Naive Added ###############

sc <- merge(d.sc[[1]], y = unlist(d.sc)[2:length(d.sc)]) # 14894 20000
sc@meta.data$orig.ident <- factor(sc@meta.data$orig.ident, levels = c("chromium033", "chromium034", "chromium035", "chromium040"))
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 5000)

sc <- ScaleData(sc, features = rownames(sc), vars.to.regress = c("nFeature_RNA", "nCount_RNA")) 
sc <- RunPCA(sc, features = VariableFeatures(object = sc)) 

#saveRDS(sc, paste('Data/mergedsc_runPCA_completed_with_naive_0406.rds', sep='')) 
sc = readRDS('Data/mergedsc_runPCA_completed_with_naive_0406.rds')

## Figure 
p = VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = 'orig.ident', ncol = 3, pt.size = 0.001)
#save_plot(paste('./Figures/vlnplot.QC_mergedsc_with_naive_0406.png',sep='.'), p, base_height = 7.5, base_width = 15)

plot1 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", pt.size = 0.5)
plot2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p = plot1 + plot2
#save_plot(paste('./Figures/scatterplot.QC_mergedsc_with_naive_0406.png',sep='.'), p, base_height = 6, base_width = 10)

options(repr.plot.height = 3, repr.plot.width = 5)
sc1 <- RunHarmony(sc, "orig.ident", plot_convergence = T) # 9 iterations
 
# Compare the Plot - before & after harmony

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = sc, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sc, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
#save_plot(paste('./Figures/dim.vimplot_mergedsc_with_naive_0406.png',sep='.'), last_plot(), base_height = 5, base_width = 12)
p1 <- DimPlot(object = sc1, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sc1, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
#save_plot(paste('./Figures/dim.vimplot_Harmony_with_naive_0406.png',sep='.'), last_plot(), base_height = 5, base_width = 12)

## Perplexity : determine clustering strength (higher value -> more distinct clusters)
ndims = 10
sc1 <- RunTSNE(sc1, reduction = "harmony", dims = 1:ndims, min.dist = 0.1, perplexity = 100) %>% # min.dist?
  FindNeighbors(reduction = "harmony", dims = 1:ndims) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

#saveRDS(sc1, paste('Data/Harmony_cluster_0to16_per100_with_naive_0406.rds', sep=''))
sc1 <- readRDS('Data/Harmony_cluster_0to16_per100_with_naive_0406.rds')

DimPlot(sc1, reduction = "tsne", label = TRUE, pt.size = .1)
save_plot(paste('./Figures/TSNE.Harmony_with_naive_all_0406.png',sep='.'), last_plot(), base_height = 6, base_width = 7)
DimPlot(sc1, reduction = "tsne", label = TRUE, pt.size = .1, split.by = "orig.ident")
save_plot(paste('./Figures/TSNE.Harmony_with_naive_split_0406.png',sep='.'), last_plot(), base_height =5, base_width = 15)

ndims = 10
sc2 <- RunUMAP(sc1, reduction = "harmony", dims = 1:ndims, min.dist = 0.1, perplexity = 100) %>% # min.dist?
  FindNeighbors(reduction = "harmony", dims = 1:ndims) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

saveRDS(sc1, paste('Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0630.rds', sep=''))


DimPlot(sc2, reduction = "umap", label = TRUE, pt.size = .1)
save_plot(paste('./Figures/UMAP.Harmony_with_naive_all_0629.png',sep='.'), last_plot(), base_height = 6, base_width = 7)
DimPlot(sc2, reduction = "umap", label = TRUE, pt.size = .1, split.by = "orig.ident")
save_plot(paste('./Figures/UMAP.Harmony_with_naive_split_0629.png',sep='.'), last_plot(), base_height =5, base_width = 15)



## Plot Cluster Tree
sc1 <- BuildClusterTree(sc1, reorder.numeric = T)
Tool(object = sc1, slot = 'BuildClusterTree')
PlotClusterTree(sc1)
save_plot(paste('./Figures/plot.clustertree_harmony_withnaive_0406.png',sep='.'), PlotClusterTree(sc1), base_height = 7, base_width = 7)

## Frequency 
prop.table(table(Idents(sc1)))
table(Idents(sc1), sc1$orig.ident)
d = as.data.frame(prop.table(table(Idents(sc1), sc1$orig.ident)))
write.table(d, 'Tables/table.frequecy.harmony.with_naive.0406.txt', sep='\t', quote = F, row.names = F, col.names = T)

## Marker
# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(sc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, paste('Tables/table.FindAllMarkers_with_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)

t = markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
write.table(t, paste('Tables/table.Top15FindAllMarkers_with_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)

