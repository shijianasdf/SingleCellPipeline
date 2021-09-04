# 4. Manual annotation

### Retrieving marker genes
query_seur = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")
markers_seur = read.delim("Tables/table.FindAllMarkers_with_naive.UMAP.Harmony.per100.0703.txt")

require(dplyr)
# Retrieve the top 5 marker genes per cluster
# Use whichever genes have the highest values under the AVG_LOG column
top5 <- markers_seur %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_seur), value = TRUE)),
                   n = 5)
# Create the dot plot
Seurat::DotPlot(query_seur, features = unique(top5$gene)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
  Seurat::NoLegend()
# Create the heatmap
Seurat::DoHeatmap(query_seur, features = unique(top5$gene)) +
  Seurat::NoLegend() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8)) # Error exists.


### Pathway analysis
# First get marker genes through cerebro
query_seur <- cerebroApp::getMarkerGenes(query_seur,
                                         groups = c('seurat_clusters'),
                                         assay = "RNA",
                                         organism = "mm",
                                         name = 'cerebro_seurat',
                                         only_pos = TRUE)

# Get enriched pathways through cerebro
query_seur <- cerebroApp::getEnrichedPathways(query_seur,
                                              # databases = c("GO_Biological_Process_2018",
                                              #               "GO_Cellular_Component_2018",
                                              #               "GO_Molecular_Function_2018",
                                              #               "KEGG_2016",
                                              #               "WikiPathways_2016",
                                              #               "Reactome_2016",
                                              #               "Panther_2016",
                                              #               "Mouse_Gene_Atlas"),
                                              adj_p_cutoff = 0.05,
                                              max_terms = 1000,
                                              URL_API = "http://amp.pharm.mssm.edu/Enrichr/enrich")

# Enriched pathways are stored in the following location:
query_seur@misc$enriched_pathways

databases = c("GO_Biological_Process_2018","GO_Cellular_Component_2018","GO_Molecular_Function_2018",
              "KEGG_2016", "WikiPathways_2016","Reactome_2016", "Panther_2016", "Mouse_Gene_Atlas")
try(
  data_from_enrichr <- .send_enrichr_query(
    marker_genes_current_group_level,
    databases = databases, URL_API = URL_API
  )
)



## gprofileR ?







