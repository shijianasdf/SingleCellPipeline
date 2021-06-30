options(stringsAsFactors = F)
library(ensembldb)
library(dplyr)
library(AnnotationHub)
library(tidyverse)

#########################
### Annotation Data

# Connect to AnnotationHub and access Ensembl db - Mus musculus 
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)

# Acquire the latest annotation files and download appropriate Ensembldb database
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
edb <- ah[[id]]

# Extract gene-level information from database and select annotations of interest
annotations <- genes(edb, return.type = "data.frame")
annotations <- annotations %>% 
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

##############################################
### Combining markers with gene descriptions 
#annotated_dataframe <- dataframe %>% 
#  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
################################################