library(parallel)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(readxl)
source("function_heatmap_SYK.R")

annot.f3 = read.csv("Data/annot_Fig_3.csv") # 17460 obs.
annot.f5 = read.csv("Data/annot_Fig_S5.csv") # 19510 obs.
harmony = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0630.rds")
data.f3 = readRDS("Data/reference_f3_0324.rds")
data.f5 = readRDS("Data/reference_fs5_0324.rds")


prop.table(table(Idents(data.f3)))
d.annot = list()
for (i in 0:10){
  d = as.data.frame(WhichCells(data.f3, idents = i))
  colnames(d) = "cell"
  cluster = paste("cluster", i, sep = ".")
  d.annot[[cluster]] = merge(d, annot.f3, by = "cell")
}

df3 = read.delim("Tables/Reference_Markers/table.FindAllMarkers.PC10.data.f3.res0.5.txt")
df3 = df3[,1:7]
df3$cluster = ifelse(df3$cluster == "0", "Non-migratory cDC2",
                     ifelse(df3$cluster == "1", "Migratory cDC2",
                            ifelse(df3$cluster == "2", "MC", 
                                   ifelse(df3$cluster == "3", "MC",
                                          ifelse(df3$cluster == "4", "Migratory cDC1",
                                                 ifelse(df3$cluster == "5", "Non-migratory cDC1",
                                                        ifelse(df3$cluster == "6", "inf-cDC2",
                                                               ifelse(df3$cluster == "7", "Non-migratory cDC2, pDC",
                                                                      ifelse(df3$cluster == "8", "Proliferating DC",
                                                                             ifelse(df3$cluster == "10", "MC",
                                                                                    ifelse(df3$cluster == "11", "Migratory cDC2", NA)))))))))))
df3 = df3[!is.na(df3$cluster),]
df3 = df3[,6:7]
colnames(df3) = c("cluster_ref", "gene")
df3 = df3 %>% unique() # 5620 -> 5037

harmony.marker = read.delim(paste('Tables/MarkerGenes_Harmony/table.FindAllMarkers_with_naive.UMAP.Harmony.per100.0630.txt', sep='.'))
df3.2 = df3 %>% group_by(gene) %>%
  mutate(cluster_ref, paste(cluster_ref, collapse = ", "))
df3.2 = df3.2[,2:3] 
colnames(df3.2) = c("gene", "cluster_ref")
df3.2 = df3.2 %>% unique

s <- strsplit(df3.2$cluster_ref, split = ", ")
df3.2 = data.frame(gene = rep(df3.2$gene, sapply(s, length)), cluster_ref = unlist(s))
df3.2 = df3.2 %>% unique() 
df3.2 = df3.2 %>% group_by(gene) %>%
  mutate(cluster_ref, paste(cluster_ref, collapse = ", "))

df3.2 = df3.2[,c(1,3)] 
colnames(df3.2) = c("gene", "cluster_ref")
df3.2 = df3.2 %>% unique

harmony.marker.annot = merge(df3.2, harmony.marker, by = "gene")
harmony.marker.annot = harmony.marker.annot %>% arrange(p_val) %>% arrange(cluster)

write.table(harmony.marker.annot, "Tables/MarkerGenes_UMAP.harmony_withreference_f3_0701.txt", sep='\t', quote = F, row.names = F, col.names = T)



### S5
prop.table(table(Idents(data.f5)))
d.annot = list()
for (i in 0:10){
  d = as.data.frame(WhichCells(data.f5, idents = i))
  colnames(d) = "cell"
  cluster = paste("cluster", i, sep = ".")
  d.annot[[cluster]] = merge(d, annot.f5, by = "cell")
}

df5 = read.delim("Tables/Reference_Markers/table.FindAllMarkers.PC10.data.fs5.res0.5.txt")
df5 = df5[,1:7]
df5$cluster = ifelse(df5$cluster == "0", "MC",
                     ifelse(df5$cluster == "1", "MC",
                            ifelse(df5$cluster == "2", "Non-migratory cDC2", 
                                   ifelse(df5$cluster == "3", "Migratory cDC2",
                                          ifelse(df5$cluster == "4", "Non-migratory cDC1",
                                                 ifelse(df5$cluster == "5", "Migratory cDC1",
                                                        ifelse(df5$cluster == "6", "inf-cDC2, Migratory cDC2",
                                                               ifelse(df5$cluster == "7", "Non-migratory cDC2",
                                                                      ifelse(df5$cluster == "8", "Proliferating DC",
                                                                             ifelse(df5$cluster == "9", "Proliferating DC",NA))))))))))
df5 = df5[!is.na(df5$cluster),] # 5699 -> 5604
df5 = df5[,6:7]
colnames(df5) = c("cluster_ref", "gene")
df5 = df5 %>% unique() 


df5.2 = df5 %>% group_by(gene) %>%
  mutate(cluster_ref, paste(cluster_ref, collapse = ", "))
df5.2 = df5.2[,2:3] 
colnames(df5.2) = c("gene", "cluster_ref")
df5.2 = df5.2 %>% unique

s <- strsplit(df5.2$cluster_ref, split = ", ")
df5.2 = data.frame(gene = rep(df5.2$gene, sapply(s, length)), cluster_ref = unlist(s))
df5.2 = df5.2 %>% unique() 
df5.2 = df5.2 %>% group_by(gene) %>%
  mutate(cluster_ref, paste(cluster_ref, collapse = ", "))

df5.2 = df5.2[,c(1,3)] 
colnames(df5.2) = c("gene", "cluster_ref")
df5.2 = df5.2 %>% unique
harmony.marker.annot = merge(df5.2, harmony.marker, by = "gene")
harmony.marker.annot = harmony.marker.annot %>% arrange(p_val) %>% arrange(cluster)

write.table(harmony.marker.annot, "Tables/MarkerGenes_UMAP.harmony_withreference_fs5_0701.txt", sep='\t', quote = F, row.names = F, col.names = T)

df = df %>% group_by(gene) %>% mutate(cluster_ref, paste(cluster_ref, collapse = ", "))
df = df[,2:3] 
colnames(df) = c("gene", "cluster_ref")
df = df %>% unique


harmony.marker.annot = merge(df, harmony.marker, by = "gene")
harmony.marker.annot = harmony.marker.annot %>% arrange(p_val) %>% arrange(cluster)

write.table(harmony.marker.annot, "Tables/MarkerGenes_UMAP.harmony_withreference_merged_0701.txt", sep='\t', quote = F, row.names = F, col.names = T)


### Fisher Test

df3 = read.delim("Tables/Reference_Markers/table.FindAllMarkers.PC10.data.f3.res0.5.txt")
df3 = df3[,1:7]
df3$cluster = ifelse(df3$cluster == "0", "Non-migratory cDC2",
                     ifelse(df3$cluster == "1", "Migratory cDC2",
                            ifelse(df3$cluster == "2", "MC", 
                                   ifelse(df3$cluster == "3", "MC",
                                          ifelse(df3$cluster == "4", "Migratory cDC1",
                                                 ifelse(df3$cluster == "5", "Non-migratory cDC1",
                                                        ifelse(df3$cluster == "6", "inf-cDC2",
                                                               ifelse(df3$cluster == "7", "Non-migratory cDC2, pDC",
                                                                      ifelse(df3$cluster == "8", "Proliferating DC",
                                                                             ifelse(df3$cluster == "10", "MC",
                                                                                    ifelse(df3$cluster == "11", "Migratory cDC2", NA)))))))))))
df3 = df3[!is.na(df3$cluster),]
df3 = df3[,6:7]
colnames(df3) = c("cluster_ref", "gene")
df3 = df3 %>% unique() # 5620 -> 5037

df5 = read.delim("Tables/Reference_Markers/table.FindAllMarkers.PC10.data.fs5.res0.5.txt")
df5 = df5[,1:7]
df5$cluster = ifelse(df5$cluster == "0", "MC",
                     ifelse(df5$cluster == "1", "MC",
                            ifelse(df5$cluster == "2", "Non-migratory cDC2", 
                                   ifelse(df5$cluster == "3", "Migratory cDC2",
                                          ifelse(df5$cluster == "4", "Non-migratory cDC1",
                                                 ifelse(df5$cluster == "5", "Migratory cDC1",
                                                        ifelse(df5$cluster == "6", "inf-cDC2, Migratory cDC2",
                                                               ifelse(df5$cluster == "7", "Non-migratory cDC2",
                                                                      ifelse(df5$cluster == "8", "Proliferating DC",
                                                                             ifelse(df5$cluster == "9", "Proliferating DC",NA))))))))))
df5 = df5[!is.na(df5$cluster),] # 5699 -> 5604
df5 = df5[,6:7]
colnames(df5) = c("cluster_ref", "gene")
df5 = df5 %>% unique() 

df = rbind.data.frame(df3, df5) # 9762
df = df %>% unique() # 7689
celltypes = levels(df$cluster_ref)

df.edit = df %>% filter(cluster_ref %in% c("inf-cDC2, Migratory cDC2", "Non-migratory cDC2, pDC")) # 798 obs

s <- strsplit(df.edit$cluster_ref, split = ", ")
df.edit = data.frame(gene = rep(df.edit$gene, sapply(s, length)), cluster_ref = unlist(s))

df = rbind.data.frame(df, df.edit)
df = df %>% unique() # 8667
df = df %>% filter(!cluster_ref %in% c("inf-cDC2, Migratory cDC2", "Non-migratory cDC2, pDC"))

str = df$cluster_ref %>% unique()
df.list = list()
for (i in str){
  gg = df %>% filter(df$cluster_ref == i) %>% dplyr::select(gene) %>% unique()
  df.list[[i]] <- gg %>% unlist()
}
df.list[['g_bg']] <- unique(unlist(df.list))

d = read.delim(paste('Tables/MarkerGenes_Harmony/table.FindAllMarkers_with_naive.UMAP.Harmony.per100.0630.txt', sep='.'))
## Make a list for the single cell data
sc_list = list()
for (i in 0:16){
  gg = d %>% filter(cluster == i) %>% dplyr::select(gene) %>% unique() 
  gg = gg[!is.na(gg$gene), ]
  sc_list[[paste("cluster", i, sep = "_")]] <- gg 
}
sc_list[['g_bg']] <- unique(unlist(sc_list))


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

df.fisher = run_geneEnrich(sc_list, df.list)
df.fisher = df.fisher %>% filter(df.fisher$fisher_p <= 0.05 & fisher_OR >= 1)
df.fisher$setA <- factor(df.fisher$setA, levels = c("cluster_0", "cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5",
                                                    "cluster_6", "cluster_7", "cluster_8", "cluster_9", "cluster_10", "cluster_11",
                                                    "cluster_12", "cluster_13", "cluster_14", "cluster_15", "cluster_16"))
df.fisher = df.fisher %>% arrange(setA)

write.table(df.fisher, 'Tables/UMAP.harmony_fishertest_reference_0701.txt', sep='\t', quote = F, row.names = F, col.names = T)

out <- plotHeatmap(gsA = sc_list, ## Vertical
                   gsB = df.list,
                   isver = T,
                   ishor = F,
                   orderA = names(sc_list),
                   orderB = names(df.list),
                   title = 'Bosteels and Single Cell GSEA')

####################


df.lin = read_xlsx("Data/Lin.et.al_2021_markergenes.xlsx")
df.lin = df.lin[,6:7]
df.lin$cluster_ref = ifelse(df.lin$cluster == "0", "HSC/MPP",
                             ifelse(df.lin$cluster == "1", "lymphoid progenitor",
                                    ifelse(df.lin$cluster == "2", "cDC",
                                           ifelse(df.lin$cluster == "3", "myeloid progenitor",
                                                  ifelse(df.lin$cluster == "4", "monocyte",
                                                         ifelse(df.lin$cluster == "5", "pDC", "neutrophil"))))))
df.lin = df.lin[,c(3,2)]
df.lin = df.lin %>% unique() 

str = df.lin$cluster_ref %>% unique()
df.list = list()
for (i in str){
  gg = df.lin %>% filter(df.lin$cluster_ref == i) %>% dplyr::select(gene) %>% unique()
  df.list[[i]] <- gg %>% unlist()
}
df.list[['g_bg']] <- unique(unlist(df.list))


df.fisher = run_geneEnrich(sc_list, df.list)
df.fisher = df.fisher %>% filter(df.fisher$fisher_p <= 0.05 & fisher_OR >= 1)
#df.fisher$setA <- factor(df.fisher$setA, levels = c("cluster_0", "cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5",
#                                                    "cluster_6", "cluster_7", "cluster_8", "cluster_9", "cluster_10", "cluster_11",
#                                                    "cluster_12", "cluster_13", "cluster_14", "cluster_15"))
#df.fisher = df.fisher %>% arrange(setA)

#write.table(df.fisher, 'Tables/harmony_fishertest_reference_0601.txt', sep='\t', quote = F, row.names = F, col.names = T)

out <- plotHeatmap(gsA = sc_list, ## Vertical
                   gsB = df.list,
                   isver = T,
                   ishor = F,
                   orderA = names(sc_list),
                   orderB = names(df.list),
                   title = 'Lin and Single Cell GSEA')


  






harmony.marker = read.delim(paste('Tables/MarkerGenes_Harmony/table.FindAllMarkers', 'HARMONY', 'per80', '0315', 'txt', sep='.'))
df = df %>% group_by(gene) %>%
  mutate(cluster_ref, paste(cluster_ref, collapse = ", "))
df = df[,2:3] 
colnames(df) = c("gene", "cluster_ref")
df = df %>% unique

harmony.marker.annot = merge(df, harmony.marker, by = "gene")
harmony.marker.annot = harmony.marker.annot %>% arrange(p_val) %>% arrange(cluster)

write.table(harmony.marker.annot, "Tables/MarkerGenes_harmony_withreference_mixed_0330.txt", sep='\t', quote = F, row.names = F, col.names = T)

harmony.list = list()
for (i in 0:14){
  listname = paste("cluster", i, sep = ".")
  harmony.list[[listname]] = harmony.marker.annot %>% filter(cluster == i)
}

harmony.marker.annot2 = harmony.marker.annot %>% group_by(gene, cluster) %>% 
  mutate(cluster_ref = paste(cluster_ref, collapse=", ")) 

d1 = FindMarkers(harmony, ident.1 = "4", ident.2 = "13", test.use = "MAST")
d2 = FindMarkers(harmony, ident.1 = "2", ident.2 = "11", test.use = "MAST")
