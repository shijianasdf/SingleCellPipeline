library(monocle3)
library(dplyr)
library(tidyverse)
library(Seurat)
library(cowplot)
library(SingleCellExperiment)
source("function_monocle.R") # Tomonocle3, get_earlier_principlal_node


## Harmony - with naive. UMAP
harmonyobj = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")

cds = ToMonocle3(harmonyobj, scale_all = TRUE, assay = "RNA", reduction_for_projection = "pca", UMAP_cluster_slot = NULL)

cds <- learn_graph(cds)
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

p <- p + labs(color_cells_by = "Harmony clusters")

p

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

## Harmony - with naive. UMAP Euclidean
harmonyobj = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")


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
                     root_pr_nodes=get_earliest_principal_node(cds, seurat_cluster = 5))
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

########################


##################### For Each Samples ########################

harmony_less <- subset(x = harmonyobj, subset = orig.ident == "chromium033")
harmony_cog <- subset(x = harmonyobj, subset = orig.ident == "chromium034")
harmony_sep <- subset(x = harmonyobj, subset = orig.ident == "chromium035")
harmony_naive <- subset(x = harmonyobj, subset = orig.ident == "chromium040")
samplelist = c(harmony_less, harmony_cog, harmony_sep, harmony_naive)
samplename = c("harmony_less", "harmony_cog", "harmony_sep", "harmony_naive")

interested_lineage_cds <- cds[rowData(cds)$gene_short_name %in% interested_genes, ]
interested_lineage_cds <- subset(x = interested_lineage_cds, subset = orig.ident == "chromium033")

i = 1
for (sample in samplelist){
  cds = ToMonocle3(sample, scale_all = TRUE, assay = "RNA", reduction_for_projection = "pca", UMAP_cluster_slot = NULL)
  cds <- learn_graph(cds)
  
  cds <- order_cells(cds, reduction_method = "UMAP", 
                     root_pr_nodes=get_earliest_principal_node(cds, seurat_cluster = 5))
  plot_cells(
    cds = cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )
  
  
  interested_genes <- c("Itgax", "H2-Aa", "H2-Ab1", "H2-Ea", "H2-Eb1", "Cd74", "H2-Dmb1", "Itgam")
  interested_lineage_cds <- cds[rowData(cds)$gene_short_name %in% interested_genes, ]

  p <- plot_genes_in_pseudotime(interested_lineage_cds[,colData(cds)$orig.ident == "chromium033"],
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








