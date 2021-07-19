library(monocle3)
library(dplyr)
library(tidyverse)
library(Seurat)
library(cowplot)
library(SingleCellExperiment)


## Harmony - with naive. UMAP
harmonyobj = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")

ToMonocle3 <- function(seurat_object,
                       scale_all = FALSE,
                       assay = "SCT",
                       reduction_for_projection = "PCA",
                       UMAP_cluster_slot = NULL){
  
  if(scale_all){
    message("Getting residuals for all Seurat genes in chosen assay slot and placing in scale.data")
    seurat_genes <- rownames(seurat_object[[assay]])
    remaining_genes <- setdiff(seurat_genes, rownames(seurat_object[[assay]]@scale.data))
    if(assay == "SCT"){
      seurat_object <- Seurat::GetResidual(seurat_object, features = remaining_genes, assay = assay, umi.assay = "RNA")
    } else {
      seurat_object <- Seurat::ScaleData(seurat_object, features = rownames(seurat_object[[assay]]))
    }
  }
  
  #We prep the seurat object by creating gene loadings for ALL genes in the Seurat scale.data slot. This is done to allow downstream monocle3 functions on gene_modules to work appropriately.
  message("Projecting gene loadings for all Seurat genes in scale.data slot")
  seurat_object <- Seurat::ProjectDim(seurat_object, reduction = reduction_for_projection, assay = assay)
  
  ##################
  
  message("Initializing CDS object")
  
  #Extract Seurat's log-transformed values
  expression_matrix <- Seurat::GetAssayData(seurat_object, assay = assay, slot = "counts")
  #Extract Seurat meta_data
  meta_data <- seurat_object@meta.data
  #Extract gene names from Seurat object SCT slot to make CDS
  seurat_genes <- data.frame(gene_short_name = rownames(seurat_object[[assay]]),
                             row.names = rownames(seurat_object[[assay]]))
  new_cds <- monocle3::new_cell_data_set(expression_data = expression_matrix, cell_metadata = meta_data, gene_metadata = seurat_genes)
  
  ##################
  
  message("Making an SCE object from the Seurat object to facilitate transfer of information from SCE to CDS")
  sce <- as.SingleCellExperiment(seurat_object, assay = assay)
  message("Loading in all Seurat reductions (PCA, HARMONY, UMAP, etc.) into CDS")
  SingleCellExperiment::reducedDims(new_cds) <- SingleCellExperiment::reducedDims(sce)
  message("Loading in specified Seurat assay into CDS")
  SummarizedExperiment::assays(new_cds) <- SummarizedExperiment::assays(sce)
  message("Loading in Seurat gene names into CDS")
  SummarizedExperiment::rowData(new_cds) <- SummarizedExperiment::rowData(sce)
  SummarizedExperiment::rowData(new_cds)$gene_short_name <-  row.names(new_cds)
  message("Loading in Seurat gene loadings into CDS")
  new_cds@preprocess_aux$gene_loadings <- seurat_object@reductions[[reduction_for_projection]]@feature.loadings.projected
  
  ##################
  
  message("Get user specified selected clusters (or active idents) from Seurat and load into CDS")
  if(is.null(UMAP_cluster_slot)){
    list_cluster <- Idents(seurat_object)
  } else {
    Idents(seurat_object) <- UMAP_cluster_slot
    list_cluster <- Idents(seurat_object)
  }
  new_cds@clusters[["UMAP"]]$clusters <- list_cluster
  #The next two commands are run in order to allow "order_cells" to be run in monocle3
  rownames(new_cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
  colnames(SingleCellExperiment::reducedDims(new_cds)[["UMAP"]]) <- NULL
  
  ##################
  
  message("Setting all cells as belonging to one partition (multiple partitions not supported yet)")
  recreate_partition <- c(rep(1, length(new_cds@colData@rownames)))
  names(recreate_partition) <- new_cds@colData@rownames
  recreate_partition <- as.factor(recreate_partition)
  new_cds@clusters[["UMAP"]]$partitions <- recreate_partition
  
  ##################
  message("Done")
  new_cds
}

cds = ToMonocle3(harmonyobj, scale_all = TRUE, assay = "RNA", reduction_for_projection = "pca", UMAP_cluster_slot = NULL)
cds <- learn_graph(cds)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, seurat_cluster){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == seurat_cluster)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes=get_earliest_principal_node(cds, seurat_cluster = 5))

plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)


interested_genes <- c("Itgax", "H2-Aa", "H2-Ab1", "H2-Ea", "H2-Eb1", "Cd74", "H2-Dmb1", "Itgam")
interested_lineage_cds <- cds[rowData(cds)$gene_short_name %in% interested_genes,]

p <- plot_genes_in_pseudotime(interested_lineage_cds,
                         color_cells_by="seurat_clusters",
                         cell_size = 0.75,
                         label_by_short_name = TRUE,
                         trend_formula = "~ splines::ns(pseudotime, df=3)",
                         vertical_jitter = NULL,
                         horizontal_jitter = NULL,
                         min_expr=0.5)
save_plot('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression_mhcii.cd11_cluster.png', p, base_height = 8, base_width = 10)


p <- plot_genes_in_pseudotime(interested_lineage_cds,
                         color_cells_by="pseudotime",
                         cell_size = 0.75,
                         label_by_short_name = TRUE,
                         trend_formula = "~ splines::ns(pseudotime, df=3)",
                         vertical_jitter = NULL,
                         horizontal_jitter = NULL,
                         min_expr=0.5)
save_plot('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression_mhcii.cd11_pseudotime.png', p, base_height = 8, base_width = 10)



interested_genes <- c("Flt3", "Cd80", "Cd86", "Fcgr1", "Adgre1")
interested_lineage_cds <- cds[rowData(cds)$gene_short_name %in% interested_genes,]

p <- plot_genes_in_pseudotime(interested_lineage_cds,
                         color_cells_by="seurat_clusters",
                         cell_size = 0.75,
                         label_by_short_name = TRUE,
                         trend_formula = "~ splines::ns(pseudotime, df=3)",
                         vertical_jitter = NULL,
                         horizontal_jitter = NULL,
                         min_expr=0.5)
save_plot('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression_others_cluster.png', p, base_height = 8, base_width = 10)

p <- plot_genes_in_pseudotime(interested_lineage_cds,
                         color_cells_by="pseudotime",
                         cell_size = 0.75,
                         label_by_short_name = TRUE,
                         trend_formula = "~ splines::ns(pseudotime, df=3)",
                         vertical_jitter = NULL,
                         horizontal_jitter = NULL,
                         min_expr=0.5)
save_plot('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression_others_pseudotime.png', p, base_height = 8, base_width = 10)

p <- plot_genes_in_pseudotime(interested_lineage_cds,
                         color_cells_by="seurat_clusters",
                         cell_size = 0.75,
                         label_by_short_name = TRUE,
                         trend_formula = "~ splines::ns(pseudotime, df=3)",
                         vertical_jitter = NULL,
                         horizontal_jitter = NULL,
                         min_expr=NULL)
save_plot('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression_others_cluster_nocut.png', p, base_height = 8, base_width = 10)

p <- plot_genes_in_pseudotime(interested_lineage_cds,
                         color_cells_by="pseudotime",
                         cell_size = 0.75,
                         label_by_short_name = TRUE,
                         trend_formula = "~ splines::ns(pseudotime, df=3)",
                         vertical_jitter = NULL,
                         horizontal_jitter = NULL,
                         min_expr=NULL)
save_plot('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression_others_pseudotime_nocut.png', p, base_height = 8, base_width = 10)



##################### For Each Samples ########################

harmony_less <- subset(x = harmonyobj, subset = orig.ident == "chromium033")
harmony_cog <- subset(x = harmonyobj, subset = orig.ident == "chromium034")
harmony_sep <- subset(x = harmonyobj, subset = orig.ident == "chromium035")
harmony_naive <- subset(x = harmonyobj, subset = orig.ident == "chromium040")
samplelist = c(harmony_less, harmony_cog, harmony_sep, harmony_naive)
samplename = c("harmony_less", "harmony_cog", "harmony_sep", "harmony_naive")

i = 1
for (sample in samplelist){
  cds = ToMonocle3(sample, scale_all = TRUE, assay = "RNA", reduction_for_projection = "pca", UMAP_cluster_slot = NULL)
  cds <- learn_graph(cds)
  
  cds <- order_cells(cds, reduction_method = "UMAP", 
                     root_pr_nodes=get_earliest_principal_node(cds, seurat_cluster = ifelse(i == 1, 4,
                                                                                            ifelse(i == 2, 0, 
                                                                                                   ifelse(i == 3, 6,
                                                                                                          ifelse(i == 4, 4, NA))))))
  plot_cells(
    cds = cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )
  
  
  interested_genes <- c("Itgax", "H2-Aa", "H2-Ab1", "H2-Ea", "H2-Eb1", "Cd74", "H2-Dmb1", "Itgam")
  interested_lineage_cds <- cds[rowData(cds)$gene_short_name %in% interested_genes,]
  
  p <- plot_genes_in_pseudotime(interested_lineage_cds,
                                color_cells_by="seurat_clusters",
                                cell_size = 0.75,
                                label_by_short_name = TRUE,
                                trend_formula = "~ splines::ns(pseudotime, df=3)",
                                vertical_jitter = NULL,
                                horizontal_jitter = NULL,
                                min_expr=0.5)
  save_plot(paste('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime.expression', samplename[i], 'mhcii.cd11_cluster.png', sep = "_"), p, base_height = 8, base_width = 10)
  
  
  p <- plot_genes_in_pseudotime(interested_lineage_cds,
                                color_cells_by="pseudotime",
                                cell_size = 0.75,
                                label_by_short_name = TRUE,
                                trend_formula = "~ splines::ns(pseudotime, df=3)",
                                vertical_jitter = NULL,
                                horizontal_jitter = NULL,
                                min_expr=0.5)
  save_plot(paste('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression', samplename[i], 'mhcii.cd11_pseudotime.png', sep = "_"), p, base_height = 8, base_width = 10)
  
  
  interested_genes <- c("Flt3", "Cd80", "Cd86", "Fcgr1", "Adgre1")
  interested_lineage_cds <- cds[rowData(cds)$gene_short_name %in% interested_genes,]
  
  p <- plot_genes_in_pseudotime(interested_lineage_cds,
                                color_cells_by="seurat_clusters",
                                cell_size = 0.75,
                                label_by_short_name = TRUE,
                                trend_formula = "~ splines::ns(pseudotime, df=3)",
                                vertical_jitter = NULL,
                                horizontal_jitter = NULL,
                                min_expr=0.5)
  save_plot(paste('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression', samplename[i], 'others_cluster.png', sep = "_"), p, base_height = 8, base_width = 10)
  
  p <- plot_genes_in_pseudotime(interested_lineage_cds,
                                color_cells_by="pseudotime",
                                cell_size = 0.75,
                                label_by_short_name = TRUE,
                                trend_formula = "~ splines::ns(pseudotime, df=3)",
                                vertical_jitter = NULL,
                                horizontal_jitter = NULL,
                                min_expr=0.5)
  save_plot(paste('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression', samplename[i], 'others_pseudotime.png', sep = '_'), p, base_height = 8, base_width = 10)
  
  p <- plot_genes_in_pseudotime(interested_lineage_cds,
                                color_cells_by="seurat_clusters",
                                cell_size = 0.75,
                                label_by_short_name = TRUE,
                                trend_formula = "~ splines::ns(pseudotime, df=3)",
                                vertical_jitter = NULL,
                                horizontal_jitter = NULL,
                                min_expr=NULL)
  save_plot(paste('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression', samplename[i], 'others_cluster_nocut.png', sep = "_"), p, base_height = 8, base_width = 10)
  
  p <- plot_genes_in_pseudotime(interested_lineage_cds,
                                color_cells_by="pseudotime",
                                cell_size = 0.75,
                                label_by_short_name = TRUE,
                                trend_formula = "~ splines::ns(pseudotime, df=3)",
                                vertical_jitter = NULL,
                                horizontal_jitter = NULL,
                                min_expr=NULL)
  save_plot(paste('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/pseudotime_expression', samplename[i], 'others_pseudotime_nocut.png', sep = "_"), p, base_height = 8, base_width = 10)
  i = i + 1
}
