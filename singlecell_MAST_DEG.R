library(Seurat)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

seuratobj = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0630.rds")
d = FindMarkers(seuratobj, ident.1 = "2", ident.2 = "13", test.use = "MAST")
d$gene_id = rownames(d)
d$change <- ifelse(d$p_val_adj <= 0.05, ifelse(d$avg_log2FC > 0, 'Up', 'Down'), 'None')

## volcano plot

name = "cluster 2 vs 13"
label_genes = d %>% filter(p_val_adj <= 0.05) %>% pull(gene_id)
x_range = max(abs(d$avg_log2FC)) * 1.05

p <- ggplot(d, aes(x = avg_log2FC, y = -log10(p_val_adj) )) + 
  geom_point(aes(fill = change), 
             shape = 21, alpha = 1,
             na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
  theme_bw(base_size = 9) + # change theme
  labs(title = paste('DEG analysis:', name, sep = " ")) + # Add a title
  xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # y-axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add cutoffs
  geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
  ggplot2::xlim(-x_range, x_range) + 
  geom_text_repel(d = subset(d, d$gene_id %in% label_genes),
                  aes(label = gene_id),
                  size = 3.5) + 
  scale_fill_manual(values = c("Up" = "#E64B35", 
                               "Down" = "#3182bd", 
                               "None" = "grey"))
p
