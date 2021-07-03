
library(fgsea)


## Harmony data - naive included

harmony = readRDS("Data/Harmony_cluster_0to16_per100_with_naive_0406.rds")

## Bosteels Data - df3, dfs5

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
df3 = df3[!is.na(df3$cluster),] #5620

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


## NAIVE INCLUDED
d = read.delim(paste('Tables/MarkerGenes_Harmony/table.FindAllMarkers_with_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'))
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
