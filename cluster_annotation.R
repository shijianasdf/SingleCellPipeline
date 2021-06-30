options(stringsAsFactors = F)
library(ensembldb)
library(dplyr)
library(AnnotationHub)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
setwd("~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/")


######################### Gene Function Annotation #########################

### Load Single Cell cluster marker table
helpless.markers <- read.delim(paste("Tables/MarkerGenes_Seurat/table.Top15FindAllMarkers", "PC10", "helpless", "txt", sep='.'))
coghelp.markers <- read.delim(paste("Tables/MarkerGenes_Seurat/table.Top15FindAllMarkers", "PC10", "coghelp", "txt", sep='.'))
sephelp.markers <- read.delim(paste('Tables/MarkerGenes_Seurat/table.Top15FindAllMarkers', 'PC10', 'sephelp', 'txt', sep='.'))
naive.markers <- read.delim(paste('Tables/MarkerGenes_Seurat/table.Top15FindAllMarkers', 'PC10', 'naive', 'txt', sep='.'))


### Gene Function Data - Ensembl db
# Connect to AnnotationHub and access Ensembl db - Mus musculus 
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)

# Acquire the latest annotation files and download appropriate Ensembldb database
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
edb <- ah[[id]]

# Extract gene-level information from database and select annotations of interest
annotations <- genes(edb, return.type = "data.frame")
annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)


### Combine markers with gene descriptions 
annotated_helpless.markers <- helpless.markers %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
write.table(annotated_helpless.markers, paste('Tables/table.Annotated_Top15FindAllMarkers', 'PC10', 'helpless', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)

annotated_coghelp.markers <- coghelp.markers %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
write.table(annotated_coghelp.markers, paste('Tables/table.Annotated_Top15FindAllMarkers', 'PC10', 'coghelp', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)

annotated_sephelp.markers <- sephelp.markers %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),by = c("gene" = "gene_name"))
write.table(annotated_sephelp.markers, paste('Tables/table.Annotated_Top15FindAllMarkers', 'PC10', 'sephelp', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)

annotated_naive.markers <- naive.markers %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),by = c("gene" = "gene_name"))
write.table(annotated_naive.markers, paste('Tables/MarkerGenes_Seurat/table.Annotated_Top15FindAllMarkers', 'PC10', 'naive', 'txt', sep='.'), 
            sep='\t', quote = F, row.names = F, col.names = T)




########################### cellmarker & panglaodb ###########################

## Made cellmarker column with google spreadsheet
cellmarker = read_tsv("Resources/Single Cell Annotation - Annotated_Top15Markers.sephelp.tsv")

## Panglaodb was downloaed from website - filter immune system, Mus musculus
panglaodb = read_tsv("Resources/PanglaoDB_markers_27_Mar_2020.tsv", col_names = T)
pldb = panglaodb %>% dplyr::filter(organ == "Immune system")
pldb = pldb %>% dplyr::filter(!species == "Hs")
pldb = pldb[,2:3]
colnames(pldb) = c("gene", "panglaoDB")
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
pldb$gene <- tolower(pldb$gene)
pldb$gene <- firstup(pldb$gene)

pldb.merged <- pldb %>% dplyr::group_by(gene) %>%
  dplyr::summarise(panglaoDB = paste(panglaoDB, collapse = ", "))

## Merge marker table with pldb
d = read.delim('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Tables/MarkerGenes_Seurat/table.Annotated_Top15FindAllMarkers.PC10.sephelp.txt')
d = d %>% left_join(y = unique(pldb.merged[,c("gene", "panglaoDB")]), by = c("gene" = "gene"))

## Merge table with cellmarker
d = cbind.data.frame(d, cellmarker$CellMarker)
d <- rename(d, "cellmarker" = "cellmarker$CellMarker")

write.table(d, "Tables/table.Annotated_Top15FindAllMarkers_CellMarker_PanglaoDB.PC10.sephelp.txt", 
            sep='\t', quote = F, row.names = F, col.names = T)
