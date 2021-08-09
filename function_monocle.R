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


plot_genes_in_pseudotime_SYK <-function(cds_subset,
                                    min_expr=NULL,
                                    cell_size=0.75,
                                    nrow=NULL,
                                    ncol=1,
                                    panel_order=NULL,
                                    color_cells_by="pseudotime",
                                    trend_formula="~ splines::ns(pseudotime, df=3)",
                                    label_by_short_name=TRUE,
                                    vertical_jitter=NULL,
                                    horizontal_jitter=NULL){
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  tryCatch({pseudotime(cds_subset)}, error = function(x) {
    stop(paste("No pseudotime calculated. Must call order_cells first."))})
  colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
  if(!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  
  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }
  
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)),
                            msg = paste("When label_by_short_name = TRUE,",
                                        "rowData must have a column of gene",
                                        "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", "partition") |
                            color_cells_by %in% names(colData(cds_subset)),
                          msg = paste("color_cells_by must be a column in the",
                                      "colData table."))
  
  if(!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in%
                                    rowData(cds_subset)$gene_short_name))
    } else {
      assertthat::assert_that(all(panel_order %in%
                                    row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100,
                          msg = paste("cds_subset has more than 100 genes -",
                                      "pass only the subset of the CDS to be",
                                      "plotted."))
  
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  assertthat::assert_that("pseudotime" %in% names(colData(cds_subset)),
                          msg = paste("pseudotime must be a column in",
                                      "colData. Please run order_cells",
                                      "before running",
                                      "plot_genes_in_pseudotime."))
  if(!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  assertthat::assert_that(!is.null(size_factors(cds_subset)))
  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }
  
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)),
                            msg = paste("When label_by_short_name = TRUE,",
                                        "rowData must have a column of gene",
                                        "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", "partition") |
                            color_cells_by %in% names(colData(cds_subset)),
                          msg = paste("color_cells_by must be a column in the",
                                      "colData table."))
  
  if(!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in%
                                    rowData(cds_subset)$gene_short_name))
    } else {
      assertthat::assert_that(all(panel_order %in%
                                    row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100,
                          msg = paste("cds_subset has more than 100 genes -",
                                      "pass only the subset of the CDS to be",
                                      "plotted."))
  
  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[,is.finite(colData(cds_subset)$pseudotime)]
  
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  
  if (is.null(min_expr)) {
    min_expr <- 0
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_colData, by.x = "Cell", by.y = "row.names")
  
  cds_exprs$adjusted_expression <- cds_exprs$expression
  
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  
  new_data <- data.frame(pseudotime = colData(cds_subset)$pseudotime)
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)
  
  model_expectation <- model_predictions(model_tbl,
                                         new_data = colData(cds_subset))
  
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell),
                             function(x) {
                               data.frame(
                                 "expectation"=model_expectation[x$f_id,
                                                                 x$Cell])
                             })
  cds_exprs <- merge(cds_exprs, expectation)
  
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (!is.null(panel_order)) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  q <- ggplot(aes(pseudotime, expression), data = cds_exprs)
  
  
  if (!is.null(color_cells_by)) {
    q <- q + geom_point(aes_string(color = color_cells_by),
                        size = I(cell_size),
                        position=position_jitter(horizontal_jitter,
                                                 vertical_jitter))
    if (class(colData(cds_subset)[,color_cells_by]) == "numeric"){
      q <- q + viridis::scale_color_viridis(option="C")
    }
  }
  else {
    q <- q + geom_point(size = I(cell_size),
                        position=position_jitter(horizontal_jitter,
                                                 vertical_jitter))
  }
  
  q <- q + geom_line(aes(x = pseudotime, y = expectation), data = cds_exprs)
  
  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow,
                                        ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }
  
  q <- q + ylab("Expression")
  
  q <- q + xlab("pseudotime")
  q <- q + monocle_theme_opts()
  q
}



monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())+
    theme(legend.title = "Harmony clusters")
}




