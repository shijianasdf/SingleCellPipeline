options(stringsAsFactors = F)
library(readxl)
library(tidyverse)
library(dplyr)
library(writexl)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(fgsea)
setwd("~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/")

pp <- fgsea::gmtPathways('Resources/c8.all.v7.2.symbols.gmt')

## Single cell samples
sample = c("helpless", "coghelp", "sephelp", "naive")
for (outtag in sample){
  print(paste("==== fgsea analysis for group :", outtag, "===="))
  d = read.delim(paste("Tables/MarkerGenes_Seurat/table.FindAllMarkers.PC10", outtag, "txt", sep = "."))
  fgseaRes = list()
  for (i in 0:length(unique(d$cluster))){
    gg = d %>% dplyr::filter(cluster == i) %>% dplyr::select(gene, avg_log2FC) %>% unique()
    gg = gg[order(gg$avg_log2FC),]
    gg = gg[!is.na(gg$gene), ]
    
    gg1 = gg %>% dplyr::filter(!is.na(avg_log2FC)) %>% pull(avg_log2FC)
    names(gg1) <- gg %>% dplyr::filter(!is.na(avg_log2FC)) %>% pull(gene) %>% toupper()
    
    fgseaRes[[paste("cluster", i, sep = "_")]] = fgsea(pathways = pp, stats = gg1, scoreType = "pos", eps = 0)
    fgseaRes[[paste("cluster", i, sep = "_")]] = fgseaRes[[paste("cluster", i, sep = "_")]] %>% dplyr::filter(padj <= 0.1 & size >= 10)
    if (nrow(fgseaRes[[paste("cluster", i, sep = "_")]]) > 0){
      dd = as.data.frame(fgseaRes[[paste("cluster", i, sep = "_")]])
      dd = dd[order(dd$NES),]
      dd$leadingEdge <- vapply(dd$leadingEdge, paste, collapse = ", ", character(1L))
      dd = dd %>% arrange(padj)
      write.table(dd, paste("Tables/Seurat_FGSEA/table_fgsea_gmt.c8.all", outtag, "cluster", i, "0404.txt", sep = "_"), quote=F, sep='\t', row.names = F, col.names = T)
      dd2 = dd[grep("HAY", dd$pathway),] %>% dplyr::filter(pval <= 0.05)
      write.table(dd2, paste("Tables/Seurat_FGSEA/table_fgsea_gmt.c8.HAY.BM", outtag, "cluster", i, "0404.txt", sep = "_"), quote=F, sep='\t', row.names = F, col.names = T)
     }
  }
}

## Harmony output

# Without filtering - default

## NAIVE EXCLUDED
d = read.delim(paste('Tables/table.FindAllMarkers_without_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'))
fgseaRes = list()
for (i in 0:length(unique(d$cluster))){
  print(paste("==== fgsea analysis for cluster number :", i, "===="))
  gg = d %>% filter(cluster == i) %>% dplyr::select(gene, avg_log2FC) %>% unique()
  gg = gg[order(gg$avg_log2FC),]
  gg = gg[!is.na(gg$gene), ]
  
  gg1 = gg %>% filter(!is.na(avg_log2FC)) %>% pull(avg_log2FC)
  names(gg1) <- gg %>% filter(!is.na(avg_log2FC)) %>% pull(gene) %>% toupper()
  
  fgseaRes[[paste("cluster", i, sep = "_")]] = fgsea(pathways = pp, stats = gg1, scoreType = "pos", eps = 0)
  fgseaRes[[paste("cluster", i, sep = "_")]] = fgseaRes[[paste("cluster", i, sep = "_")]] %>% filter(padj <= 0.1 & size >= 10)
  if (nrow(fgseaRes[[paste("cluster", i, sep = "_")]]) > 0){
    dd = as.data.frame(fgseaRes[[paste("cluster", i, sep = "_")]])
    dd = dd[order(dd$NES),]
    dd$leadingEdge <- vapply(dd$leadingEdge, paste, collapse = ", ", character(1L))
    dd = dd %>% arrange(padj)
    write.table(dd, paste("Tables/table_fgsea_gmt.c8.all.HARMONY_wihtoutnaive",  "cluster", i, "0413.txt", sep = "_"), quote=F, sep='\t', row.names = F, col.names = T)
    dd2 = dd[grep("HAY", dd$pathway),] %>% dplyr::filter(pval <= 0.05)
    write.table(dd2, paste("Tables/table_fgsea_gmt.c8.all.HAY_BM.HARMONY_wihtoutnaive", "cluster", i, "0404.txt", sep = "_"), quote=F, sep='\t', row.names = F, col.names = T)
    }
}

## NAIVE INCLUDED
d = read.delim(paste('Tables/table.FindAllMarkers_with_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'))
fgseaRes = list()
for (i in 0:length(unique(d$cluster))){
  print(paste("==== fgsea analysis for cluster number :", i, "===="))
  gg = d %>% filter(cluster == i) %>% dplyr::select(gene, avg_log2FC) %>% unique()
  gg = gg[order(gg$avg_log2FC),]
  gg = gg[!is.na(gg$gene), ]
  
  gg1 = gg %>% filter(!is.na(avg_log2FC)) %>% pull(avg_log2FC)
  names(gg1) <- gg %>% filter(!is.na(avg_log2FC)) %>% pull(gene) %>% toupper()
  
  fgseaRes[[paste("cluster", i, sep = "_")]] = fgsea(pathways = pp, stats = gg1, scoreType = "pos", eps = 0)
  fgseaRes[[paste("cluster", i, sep = "_")]] = fgseaRes[[paste("cluster", i, sep = "_")]] %>% filter(padj <= 0.1 & size >= 10)
  if (nrow(fgseaRes[[paste("cluster", i, sep = "_")]]) > 0){
    dd = as.data.frame(fgseaRes[[paste("cluster", i, sep = "_")]])
    dd = dd[order(dd$NES),]
    dd$leadingEdge <- vapply(dd$leadingEdge, paste, collapse = ", ", character(1L))
    dd = dd %>% arrange(padj)
    write.table(dd, paste("Tables/table_fgsea_gmt.c8.all.HARMONY_default",  "cluster", i, "0413.txt", sep = "_"), quote=F, sep='\t', row.names = F, col.names = T)
    dd2 = dd[grep("HAY", dd$pathway),] %>% dplyr::filter(pval <= 0.05)
    write.table(dd2, paste("Tables/table_fgsea_gmt.c8.all.HAY_BM.HARMONY_withnaive", "cluster", i, "0404.txt", sep = "_"), quote=F, sep='\t', row.names = F, col.names = T)
    }
}

# # With no lfc threshold - X
# d = read.delim(paste('Tables/table.FindAllMarkers_no_lfcthreshold', 'per80', 'HARMONY','0315', 'txt', sep='.'))
# fgseaRes = list()
# for (i in 0:length(unique(d$cluster))){
#   print(paste("==== fgsea analysis for cluster number :", i, "===="))
#   gg = d %>% filter(cluster == i) %>% dplyr::select(gene, avg_log2FC) %>% unique()
#   gg = gg[order(gg$avg_log2FC),]
#   gg = gg[!is.na(gg$gene), ]
#   
#   gg1 = gg %>% filter(!is.na(avg_log2FC)) %>% pull(avg_log2FC)
#   names(gg1) <- gg %>% filter(!is.na(avg_log2FC)) %>% pull(gene) %>% toupper()
#   
#   fgseaRes[[paste("cluster", i, sep = "_")]] = fgsea(pathways = pp, stats = gg1, scoreType = "pos", eps = 0)
#   fgseaRes[[paste("cluster", i, sep = "_")]] = fgseaRes[[paste("cluster", i, sep = "_")]] %>% filter(padj <= 0.1 & size >= 10)
#   if (nrow(fgseaRes[[paste("cluster", i, sep = "_")]]) > 0){
#     dd = as.data.frame(fgseaRes[[paste("cluster", i, sep = "_")]])
#     dd = dd[order(dd$NES),]
#     dd$leadingEdge <- vapply(dd$leadingEdge, paste, collapse = ", ", character(1L))
#     dd = dd %>% arrange(padj)
#     write.table(dd, paste("Tables/table_fgsea_gmt.c8.all.HARMONY_nolfc",  "cluster", i, "0315.txt", sep = "_"), quote=F, sep='\t', row.names = F, col.names = T)
#   }
# }
# saveRDS(fgseaRes, file = 'Data/fgsea_Harmony_nolfc_0315.Rdata')

