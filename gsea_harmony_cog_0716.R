options(stringsAsFactors = F)
library(readxl)
library(tidyverse)
library(dplyr)
library(writexl)
setwd("~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/")
source("function_heatmap_SYK.R")

## Harmony output
d = read.delim('Tables/MarkerGenes_Harmony/table.FindAllMarkers_with_naive.UMAP.Harmony.per100.0630.txt')

## Make a list for the single cell data
sc_list = list()
for (i in 0:16){
  gg = d %>% filter(cluster == i) %>% dplyr::select(gene) %>% unique() 
  gg = gg[!is.na(gg$gene), ]
  sc_list[[paste("har", "cluster", i, sep = "_")]] <- gg %>% toupper()
}
sc_list[['g_bg']] <- unique(unlist(sc_list))


cog = read.delim('Tables/MarkerGenes_Seurat/table.FindAllMarkers.PC10.coghelp.res1.0.txt')
cog_list = list()
for (i in 0:16){
  gg = cog %>% filter(cluster == i) %>% dplyr::select(gene) %>% unique() 
  gg = gg[!is.na(gg$gene), ]
  cog_list[[paste("cog", "cluster", i, sep = "_")]] <- gg %>% toupper()
}
cog_list[['g_bg']] <- unique(unlist(cog_list))

naive = read.delim('Tables/MarkerGenes_Seurat/table.FindAllMarkers.PC10.naive.txt')
naive_list = list()
for (i in 0:10){
  gg = naive %>% filter(cluster == i) %>% dplyr::select(gene) %>% unique() 
  gg = gg[!is.na(gg$gene), ]
  naive_list[[paste("naive", "cluster", i, sep = "_")]] <- gg %>% toupper()
}
naive_list[['g_bg']] <- unique(unlist(naive_list))


# Function: Run Fisher
run_geneEnrich <- function(genesetA, genesetB){
  # genesetA: Gene list of interest (e.g. single cell cluster genes, eGenes) -> d
  # genesetB: Gene list for testing (e.g. gene set for autism neurobiology) -> pp
  # genesetA_bg: Background gene list of interest
  # genesetB_bg: Background gene list of a cell type
  
  # Find the background gene list
  genesetA_bg = sort(unique(unlist(genesetA)))
  genesetB_bg = sort(unique(unlist(genesetB)))
  
  print (length(genesetA_bg)) # 4269
  print (length(genesetB_bg)) # 12396
  
  # Run fisher exact test
  cols = c('setA', 'setB', 
           'setA_size', 'setB_size', 
           'setA_size_background', 'setB_size_background', 
           'setA_size_hits', 
           'fisher_OR', 'fisher_p', 'overlaps')
  res = as.data.frame(matrix(nrow=0, ncol=length(cols)))
  for (genesetB_1 in names(genesetB)){
    print (genesetB_1)
    genesetB1 = genesetB[[genesetB_1]]
    genesetB_bg1 = genesetB_bg[!(genesetB_bg %in% genesetB1)]
    for (genesetA_1 in names(genesetA)){
      genesetA1 = genesetA[[genesetA_1]]
      genesetA_bg1 = genesetA_bg[!(genesetA_bg %in% genesetA1)]
      
      # create a matrix
      mat = matrix( c( sum(genesetB1 %in% genesetA1), sum(genesetA_bg1 %in% genesetB1),
                       sum(genesetB_bg1 %in% genesetA1), sum(genesetA_bg1 %in% genesetB_bg1)), 
                    ncol=2   )
      
      fis = fisher.test(mat)
      
      out = as.data.frame(matrix(c(genesetA_1, genesetB_1,
                                   length(genesetA1), length(genesetB1), 
                                   length(genesetA_bg), length(genesetB_bg), 
                                   sum(genesetB1 %in% genesetA1), 
                                   fis$estimate, fis$p.value, 
                                   paste(sort(intersect(genesetB1, genesetA1)), collapse=', ' )) , ncol=length(cols)))
      colnames(out) = cols
      res = rbind.data.frame(res, out)
    }
  }
  # update the type
  res$fisher_OR <- as.numeric(res$fisher_OR)
  res$fisher_p <- as.numeric(res$fisher_p)
  return(res)
}

df = run_geneEnrich(sc_list, cog_list)
df = run_geneEnrich(sc_list, naive_list)

# ###########################
# str = df$setB %>% unique
# pp_list = list()
# for (i in str){
#   pp_list[[i]] = pp[[i]]
# }
# ###########################
# 
# df = df %>% filter(grepl("HAY", setB))
# df = df[df$setA != 'g_bg' & df$setB != 'g_bg',]
# df = df %>% filter(df$fisher_p <= 0.05 & df$fisher_OR >= 1.00)
# df = df %>% arrange(fisher_p)
# df = df %>% arrange(setA)

# write.table(df, 'Tables/harmony.without_naive_fishertest.msigdb_hay.0413.txt', sep='\t', 
#             quote = F, row.names = F, col.names = T)

out <- plotHeatmap(gsA = sc_list, 
                   gsB = cog_list,
                   isver = T,
                   ishor = F,
                   orderA = names(sc_list),
                   orderB = names(cog_list),
                   title = 'Harmony and coghelp GSEA')


out <- plotHeatmap(gsA = sc_list, 
                   gsB = naive_list,
                   isver = T,
                   ishor = F,
                   orderA = names(sc_list),
                   orderB = names(naive_list),
                   title = 'Harmony and naive GSEA')






options(stringsAsFactors = F)
library(readxl)
library(tidyverse)
library(dplyr)
library(writexl)
library(fgsea)
setwd("~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/")


## Harmony output
d = read.delim('Tables/MarkerGenes_Harmony/table.FindAllMarkers_with_naive.UMAP.Harmony.per100.0630.txt')

## Make a list for the single cell data
sc_list = list()
for (i in 0:16){
  gg = d %>% filter(cluster == i) %>% dplyr::select(gene) %>% unique() 
  gg = gg[!is.na(gg$gene), ]
  sc_list[[paste("har", "cluster", i, sep = "_")]] <- gg %>% toupper()
}
sc_list[['g_bg']] <- unique(unlist(sc_list))

cog = read.delim('Tables/MarkerGenes_Seurat/table.FindAllMarkers.PC10.coghelp.res1.0.txt')
cog_list = list()
for (i in 0:14){
  gg = cog %>% filter(cluster == i) %>% dplyr::select(gene, avg_log2FC) %>% unique()
  gg = gg[order(gg$avg_log2FC), ]
  gg = gg[!is.na(gg$gene), ]
  
  gg1 = gg %>% dplyr::filter(!is.na(avg_log2FC)) %>% pull(avg_log2FC)
  names(gg1) <- gg %>% dplyr::filter(!is.na(avg_log2FC)) %>% pull(gene) %>% toupper()
  
  cog_list[[paste("cog", "cluster", i, sep = "_")]] = fgsea(pathways = sc_list, stats = gg1, scoreType = "pos", eps = 0)
  cog_list[[paste("cog", "cluster", i, sep = "_")]] = fgseaRes[[paste("cog", "cluster", i, sep = "_")]] %>% dplyr::filter(padj <= 0.1 & size >= 10)
}


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

