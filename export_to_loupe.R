options(stringsAsFactors = F)
library(DESeq2)
library(tidyverse)
library(harmony)
library(dplyr)
library(Seurat)
setwd("~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/")

## Harmony with naive
sc1 <- readRDS('Data/Harmony_cluster_0to16_per100_with_naive_0406.rds')
tsne_sc1 = cbind("Barcode" = rownames(Embeddings(object = sc1, reduction = "tsne")), Embeddings(object = sc1, reduction = "tsne"))
#write.table(tsne_sc1, file="./Data/Harmony_cluster_0to16_per100_with_naive_0406.for_loupe.csv", sep = ",", quote = F, row.names = F, col.names = T)


## Seurat clusters
helpless <- readRDS("Data/SeuratRDS/chromium033_PC10_res0.5.rds")

# Embeddings information - UMAP Plotting
umap_helpless = cbind("Barcode" = rownames(Embeddings(object = helpless, reduction = "umap")), Embeddings(object = helpless, reduction = "umap"))
#write.table(umap_helpless, file="./Data/chromium033_PC10_res0.5.for_loupe.csv", sep = ",", quote = F, row.names = F, col.names = T)

# Category information - for clusters
dc = data.frame()
for (i in 0:12){
  dc_tmp = as.data.frame(WhichCells(helpless, idents = i))
  colnames(dc_tmp) = "barcode"
  dc_tmp$id = paste("Cluster", i, sep = " ")
  dc = rbind.data.frame(dc, dc_tmp)
}
write.csv(dc, file = "./Data/Chromium033_PC10_res0.5_cluster_annotation.csv", row.names = F)

coghelp <- readRDS("Data/SeuratRDS/chromium034_PC10_res0.5.rds")
umap_coghelp = cbind("Barcode" = rownames(Embeddings(object = coghelp, reduction = "umap")), Embeddings(object = coghelp, reduction = "umap"))
#write.table(umap_coghelp, file="./Data/chromium034_PC10_res0.5.for_loupe.csv", sep = ",", quote = F, row.names = F, col.names = T)

dc = data.frame()
for (i in 0:11){
  dc_tmp = as.data.frame(WhichCells(coghelp, idents = i))
  colnames(dc_tmp) = "barcode"
  dc_tmp$id = paste("Cluster", i, sep = " ")
  dc = rbind.data.frame(dc, dc_tmp)
}
write.csv(dc, file = "./Data/Chromium034_PC10_res0.5_cluster_annotation.csv", row.names = F)

sephelp <- readRDS("Data/SeuratRDS/chromium035_PC10_res0.5.rds")
umap_sephelp = cbind("Barcode" = rownames(Embeddings(object = sephelp, reduction = "umap")), Embeddings(object = sephelp, reduction = "umap"))
#write.table(umap_sephelp, file="./Data/chromium035_PC10_res0.5.for_loupe.csv", sep = ",", quote = F, row.names = F, col.names = T)

dc = data.frame()
for (i in 0:11){
  dc_tmp = as.data.frame(WhichCells(sephelp, idents = i))
  colnames(dc_tmp) = "barcode"
  dc_tmp$id = paste("Cluster", i, sep = " ")
  dc = rbind.data.frame(dc, dc_tmp)
}
write.csv(dc, file = "./Data/Chromium035_PC10_res0.5_cluster_annotation.csv", row.names = F)

naive <- readRDS("Data/SeuratRDS/chromium040_PC10_res0.5.rds")
umap_naive = cbind("Barcode" = rownames(Embeddings(object = naive, reduction = "umap")), Embeddings(object = naive, reduction = "umap"))
#write.table(umap_naive, file="./Data/chromium040_PC10_res0.5.for_loupe.csv", sep = ",", quote = F, row.names = F, col.names = T)

dc = data.frame()
for (i in 0:10){
  dc_tmp = as.data.frame(WhichCells(naive, idents = i))
  colnames(dc_tmp) = "barcode"
  dc_tmp$id = paste("Cluster", i, sep = " ")
  dc = rbind.data.frame(dc, dc_tmp)
}
write.csv(dc, file = "./Data/Chromium040_PC10_res0.5_cluster_annotation.csv", row.names = F)






