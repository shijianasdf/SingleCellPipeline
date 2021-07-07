library(monocle3)
library(dplyr)
library(tidyverse)
library(Seurat)
library(cowplot)
library(SingleCellExperiment)


###############################################################################################################
## Trajectory on UMAP

sample = "chromium040"
seuratobj = readRDS(paste("Data/SeuratRDS/", sample, "_PC10_res0.5.rds", sep = ""))

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

cds = ToMonocle3(seuratobj, scale_all = TRUE, assay = "RNA", reduction_for_projection = "pca", UMAP_cluster_slot = NULL)
cds <- learn_graph(cds)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, seurat_cluster="0"){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == seurat_cluster)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes=get_earliest_principal_node(cds, seurat_cluster = 4))

plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)



######################################
## Harmony - with naive. UMAP

harmonyobj = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0703.rds")

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

#cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes=get_earliest_principal_node(cds, seurat_cluster = 5))
#cds <- order_cells(cds, reduction_method = "UMAP")
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes=get_earliest_principal_node(cds, seurat_cluster = 5))
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes=c(get_earliest_principal_node(cds, seurat_cluster = 5), get_earliest_principal_node(cds, seurat_cluster = 4)))

p <- plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

save_plot('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/plot.Monocle_UMAP.harmony_cluster4&5_0706.png', p, base_height = 6, base_width = 7)


################################################################
## Harmony - divide coghelp vs naive

harmony_cog <- subset(x = harmonyobj, subset = orig.ident == "chromium034")
harmony_naive <- subset(x = harmonyobj, subset = orig.ident == "chromium040")

cds_cog = ToMonocle3(harmony_cog, scale_all = TRUE, assay = "RNA", reduction_for_projection = "pca", UMAP_cluster_slot = NULL)
cds_cog <- learn_graph(cds_cog)
cds_cog <- order_cells(cds_cog, reduction_method = "UMAP", root_pr_nodes=get_earliest_principal_node(cds_cog, seurat_cluster = 5))

p1 <- plot_cells(
  cds = cds_cog,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

cds_naive = ToMonocle3(harmony_naive, scale_all = TRUE, assay = "RNA", reduction_for_projection = "pca", UMAP_cluster_slot = NULL)
cds_naive <- learn_graph(cds_naive)
cds_naive <- order_cells(cds_naive, reduction_method = "UMAP", root_pr_nodes= get_earliest_principal_node(cds_naive, seurat_cluster = 5))

p2 <- plot_cells(
  cds = cds_naive,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

p = plot_grid(p1, p2)
save_plot('~/Dropbox/NGR_SNU_2019/scRNA_seq_SYK/Figures/plot.Monocle_UMAP.harmony_split_cluster5_0706.png', p, base_height = 6, base_width = 14)


###############################################################################################################
## Branch-like

harmonyobj = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0703.rds")

expression_matrix <- harmonyobj@assays$RNA@counts %>% as.matrix()
cell_metadata = harmonyobj@meta.data
cell_metadata <- cbind(cell_metadata, harmonyobj[["umap"]]@cell.embeddings)
gene_annotation = data.frame(gene_short_name=rownames(expression_matrix))
rownames(gene_annotation) = gene_annotation$gene_short_name

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 10) 

# Reduce dimensionality 
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Clustering
cds = cluster_cells(cds, reduction_method = "UMAP", resolution=0.7)
p = plot_cells(cds, group_cells_by="cluster", color_cells_by = 'seurat_clusters', reduction_method = 'UMAP', label_cell_groups = T)
p

##pseudotime
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes= get_earliest_principal_node(cds, seurat_cluster = 5))

p = plot_cells(cds,
               color_cells_by = "seurat_clusters",
               label_groups_by_cluster=T,
               label_leaves=T,
               label_branch_points=T, group_cells_by = 'cluster', show_trajectory_graph = T, label_cell_groups = F)

p
ggsave(paste('Figures/plot.', tag, '_pseudotime.before.', date, '.png', sep = ''), p, 
       height = 4, width = 6, limitsize = F)
cds <- order_cells(cds)
p = plot_cells(cds,
               color_cells_by = "pseudotime",
               label_cell_groups=T,
               label_leaves=T,
               label_branch_points=T,
               graph_label_size=1.5, label_groups_by_cluster=T)
p
ggsave(paste('Figures/plot.', tag, '_pseudotime.after.', date, '.png', sep = ''), p, 
       height = 4, width = 6, limitsize = F)


reducedDims(cds)$UMAP[,1] <- cell_metadata$UMAP_1
reducedDims(cds)$UMAP[,2] <- cell_metadata$UMAP_2





p<-plot_cells(cds, reduction_method = "UMAP",
              color_cells_by = "seurat_clusters", group_label_size = 3.5,
              label_groups_by_cluster = F, label_cell_groups = T)
p
ggsave(paste('Figures/plot.', tag, '_pseudotime_UMAP.before.', date, '.png', sep = ''), p, 
       height = 5, width = 6, limitsize = F)
cds <- order_cells(cds)
p<-plot_cells(cds, 
              color_cells_by = "pseudotime", group_label_size = 3.5,
              label_groups_by_cluster = FALSE)
p
ggsave(paste('Figures/plot.', tag, '_pseudotime_UMAP.after.', date, '.png', sep = ''), p, 
       height = 5, width = 7, limitsize = F)



####


cds <- learn_graph(cds)
p = plot_cells(cds,
               color_cells_by = "assigned_cell_type",
               label_groups_by_cluster=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE)
p
ggsave(paste('Figures/gRNA/plot.pseudotime.before', names(cds_list)[k], date, 'png', sep = '.'), p, 
       height = 4, width = 6, limitsize = F)
cds <- order_cells(cds)
p = plot_cells(cds,
               color_cells_by = "pseudotime",
               label_cell_groups=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE,
               graph_label_size=1.5)
p
ggsave(paste('Figures/gRNA/plot.pseudotime.after', names(cds_list)[k], date, 'png', sep = '.'), p, 
       height = 4, width = 6, limitsize = F)














