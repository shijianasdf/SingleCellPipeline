options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(Seurat)
library(cowplot)

date = "0701"

## Load Single cell Dataset
helpless = readRDS("Data/SeuratRDS/chromium033_PC10_res0.5.rds")
coghelp = readRDS("Data/SeuratRDS/chromium034_PC10_res0.5.rds")
sephelp = readRDS("Data/SeuratRDS/chromium035_PC10_res0.5.rds")
naive = readRDS("Data/SeuratRDS/chromium040_PC10_res0.5.rds")

samplelist = list(helpless, coghelp, sephelp, naive)
samplename = list("helpless", "coghelp", "sephelp", "naive")

## Load Harmony Dataset
harmony_wo = readRDS("Data/Harmony_cluster_0to16_per100_without_naive_0406.rds")
harmony_wt = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0630.rds")

samplelist = list(harmony_wo, harmony_wt)
samplename = list("harmony_wo", "harmony_wt")


###############################################
## Phenotype encoding genes
cdc1_marker <- c("Cd8a", "Xcr1", "Cd24", "Clec9a", "Ly75", "Zbtb46", "Itgae")
cdc1_marker <- str_split(cdc1_marker, ", ") %>% unlist # char type.
cdc2_marker <- c("Itgam", "Sirpa", "Cd4", "Clec10a", "Clec12a", "Zbtb46")
cdc2_marker <- str_split(cdc2_marker, ", ") %>% unlist # char type.
mo_marker <- c("Itgam", "Cd209", "Mrc1", "Mrc2", "Fcgr1", "Adgre1", "Sirpa", "Cd14")
mo_marker <- str_split(mo_marker, ", ") %>% unlist # char type.
pdc_marker <- c("Ptprc", "Bst2", "Ly49q", "Siglech", "Ccr9")
pdc_marker <- str_split(pdc_marker, ", ") %>% unlist # char type.

markerlist = list(cdc1_marker, cdc2_marker, mo_marker, pdc_marker)
markername = list("cdc1_marker", "cdc2_marker", "mo_marker", "pdc_marker")

for (i in 1:2){
  for (j in 1:3){
   # p <- VlnPlot(samplelist[[i]], markerlist[[j]], group.by = "seurat_clusters", ncol = 4)
  #  save_plot(paste("./Figures/vlnplot_", samplename[[i]], ".", markername[[j]], ".pheno_0528.pdf", sep = ""), p, base_height = 6, base_width = 14)
    p <- FeaturePlot(object = samplelist[[i]], features = markerlist[[j]], cols = c("lightgrey", "blue"), ncol = 4)
    save_plot(paste("./Figures/featureplot_", samplename[[i]], ".", markername[[j]], ".pheno_0528.pdf", sep = ""), p, base_height = 6, base_width = 14)
  }
}

j = 4
for (i in 1:2){
  p <- VlnPlot(samplelist[[i]], markerlist[[j]], group.by = "seurat_clusters", ncol = 4)
  save_plot(paste("./Figures/vlnplot_", samplename[[i]], ".", markername[[j]], ".pheno_0528.pdf", sep = ""), p, base_height = 3, base_width = 14)
  p <- FeaturePlot(object = samplelist[[i]], features = markerlist[[j]], cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_", samplename[[i]], ".", markername[[j]], ".pheno_0528.pdf", sep = ""), p, base_height = 3, base_width = 14)
}

## for UMAP harmony data
for (j in 1:3){
  p <- FeaturePlot(object = harmony_wt, features = markerlist[[j]], cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_UMAP.harmony_", markername[[j]], ".pheno_", date, ".pdf", sep = ""), p, base_height = 6, base_width = 14)
}
j = 4
p <- FeaturePlot(object = harmony_wt, features = markerlist[[j]], cols = c("lightgrey", "blue"), ncol = 4)
save_plot(paste("./Figures/featureplot_UMAP.harmony_", markername[[j]], ".pheno_", date, ".pdf", sep = ""), p, base_height = 3, base_width = 14)


###############################################
## CD11c and MHCII coding genes
marker = c("Itgax", "H2-Ea", "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "H2-Dma", "H2-DMb1")

for (i in 1:2){
  p <- VlnPlot(samplelist[[i]], marker, group.by = "seurat_clusters", ncol = 4)
  save_plot(paste("./Figures/vlnplot_", samplename[[i]], ".marker_explevel_0528.pdf", sep = ""), p, base_height = 6, base_width = 14)
  p <- FeaturePlot(object = samplelist[[i]], features = marker, cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_", samplename[[i]], ".marker_explevel_0528.pdf", sep = ""), p, base_height = 6, base_width = 14)
}

for (i in 1:2){
  p <- VlnPlot(samplelist[[i]], marker, group.by = "seurat_clusters", ncol = 4)
  save_plot(paste("./Figures/vlnplot_", samplename[[i]], ".marker_explevel_0528.pdf", sep = ""), p, base_height = 6, base_width = 14)
  p <- FeaturePlot(object = samplelist[[i]], features = marker, cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_", samplename[[i]], ".marker_explevel_0528.pdf", sep = ""), p, base_height = 6, base_width = 14)
}

p <- FeaturePlot(object = harmony_wt, features = marker, cols = c("lightgrey", "blue"), ncol = 4)
save_plot(paste("./Figures/featureplot_UMAP.harmony_marker_explevel", date, ".pdf", sep = ""), p, base_height = 6, base_width = 14)


###############################################
## Extra Marker coding genes
extra_marker <- c("Mertk", "Csf1r", "Cx3cr1", "Ccr2", "Cxcl9", "Cxcl10", "Ly6c", "Ly6g", "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2", "Il2", 
                  "Il2ra", "Il2rb", "Il2rg", "Il10ra", "Pdcd1", "Cd274a")
extra_marker <- str_split(cdc1_marker, ", ") %>% unlist # char type.

for (i in 1:4){
  p <- VlnPlot(samplelist[[i]], extra_marker, group.by = "seurat_clusters", ncol = 4)
  save_plot(paste("./Figures/vlnplot_", samplename[[i]], ".extra_marker_0522.pdf", sep = ""), p, base_height = 6, base_width = 14)
  p <- FeaturePlot(object = samplelist[[i]], features = extra_marker, cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_", samplename[[i]], ".extra_marker_0522.pdf", sep = ""), p, base_height = 6, base_width = 14)
}

for (i in 1:2){
  p <- VlnPlot(samplelist[[i]], extra_marker, group.by = "seurat_clusters", ncol = 4)
  save_plot(paste("./Figures/vlnplot_", samplename[[i]], ".extra_marker_0528.pdf", sep = ""), p, base_height = 6, base_width = 14)
  p <- FeaturePlot(object = samplelist[[i]], features = extra_marker, cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_", samplename[[i]], ".extra_marker_0528.pdf", sep = ""), p, base_height = 6, base_width = 14)
}

## for UMAP harmony data
p <- FeaturePlot(object = harmony_wt, features = extra_marker, cols = c("lightgrey", "blue"), ncol = 4)
save_plot(paste("./Figures/featureplot_UMAP.harmony_", "extra_marker_", date, ".pdf", sep = ""), p, base_height = 3, base_width = 14)



###############################################
## TLR coding genes
Tlr_marker = c("Tlr1", "Tlr2", "Tlr3", "Tlr4", "Tlr5", "Tlr6", "Tlr7", "Tlr8", "Tlr9", "Tlr11", "Tlr12", "Tlr13")

for (i in 1:4){
  p <- VlnPlot(samplelist[[i]], Tlr_marker, group.by = "seurat_clusters", ncol = 4)
  save_plot(paste("./Figures/vlnplot_", samplename[[i]], ".Tlr_marker_0522.pdf", sep = ""), p, base_height = 9, base_width = 14)
  p <- FeaturePlot(object = samplelist[[i]], features = Tlr_marker, cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_", samplename[[i]], ".Tlr_marker_0522.pdf", sep = ""), p, base_height = 9, base_width = 14)
}

for (i in 1:2){
 # p <- VlnPlot(samplelist[[i]], Tlr_marker, group.by = "seurat_clusters", ncol = 4)
#  save_plot(paste("./Figures/vlnplot_", samplename[[i]], ".", "Tlr_marker_0528.pdf", sep = ""), p, base_height = 9, base_width = 14)
  p <- FeaturePlot(object = samplelist[[i]], features = Tlr_marker, cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_", samplename[[i]], ".", "Tlr_marker_0528.pdf", sep = ""), p, base_height = 9, base_width = 14)
}

## for UMAP harmony data
p <- FeaturePlot(object = harmony_wt, features = Tlr_marker, cols = c("lightgrey", "blue"), ncol = 4)
save_plot(paste("./Figures/featureplot_UMAP.harmony_", "Tlr_marker_", date, ".pdf", sep = ""), p, base_height = 9, base_width = 14)



###############################################
## Transcription factors
cdc1_tf = c("Irf8", "Spi1", "Batf3", "Id2")
cdc2_tf = c("Irf4", "Spi1")
modc_tf = c("Irf4", "Spi1")
pdc_tf = c("Tcf4", "Zeb2")

markerlist = list(cdc1_tf, cdc2_tf, modc_tf, pdc_tf)
markername = list("cdc1_tf", "cdc2_tf", "modc_tf", "pdc_tf")

for (i in 1:4){
  for (j in 1:4){
    p <- VlnPlot(samplelist[[i]], markerlist[[j]], group.by = "seurat_clusters", ncol = 4)
    save_plot(paste("./Figures/vlnplot_", samplename[[i]], markername[[j]], "marker_0522.pdf", sep = "."), p, base_height = 3, base_width = 14)
    p <- FeaturePlot(object = samplelist[[i]], features = markerlist[[j]], cols = c("lightgrey", "blue"), ncol = 4)
    save_plot(paste("./Figures/featureplot_", samplename[[i]], markername[[j]], "marker_0522.pdf", sep = "."), p, base_height = 3, base_width = 14)
  }
}

for (i in 1:2){
  for (j in 1:4){
   # p <- VlnPlot(samplelist[[i]], markerlist[[j]], group.by = "seurat_clusters", ncol = 4)
   # save_plot(paste("./Figures/vlnplot_", samplename[[i]], markername[[j]], "marker_0528.pdf", sep = "."), p, base_height = 3, base_width = 14)
    p <- FeaturePlot(object = samplelist[[i]], features = markerlist[[j]], cols = c("lightgrey", "blue"), ncol = 4)
    save_plot(paste("./Figures/featureplot_", samplename[[i]], markername[[j]], "marker_0528.pdf", sep = "."), p, base_height = 3, base_width = 14)
  }
}

for (j in 1:4){
  p <- FeaturePlot(object = harmony_wt, features = markerlist[[j]], cols = c("lightgrey", "blue"), ncol = 4)
  save_plot(paste("./Figures/featureplot_UMAP.harmony_", markername[[j]], ".marker_", date, ".pdf", sep = ""), p, base_height = 3, base_width = 14)
}


###############################################
## Secretory Products 

protein_list = c("Il12a", "Il12b", "Il6st", "Il6ra", "Ifnar1", "Ifnar2")

for (i in 1:4){
  p <- FeaturePlot(object = samplelist[[i]], features = protein_list, cols = c("lightgrey", "blue"), ncol = 3)
  save_plot(paste("./Figures/featureplot_", samplename[[i]], ".secretory_0628.pdf", sep = ""), p, base_height = 9, base_width = 14)
}
for (i in 1:2){
    p <- FeaturePlot(object = samplelist[[i]], features = protein_list, cols = c("lightgrey", "blue"), ncol = 3)
    save_plot(paste("./Figures/featureplot_", samplename[[i]], ".secretory_0628.pdf", sep = ""), p, base_height = 9, base_width = 14)
}
