## Analyze RNA velocity quantifications stored in Seurat object using scVelo.
## https://github.com/velocyto-team/velocyto.R
options(stringsAsFactors = F)
library(dplyr)
library(tidyverse)
library(Seurat)
library(SeuratDisk) # remotes::install_github("mojaveazure/seurat-disk")
library(SeuratWrappers) # remotes::install_github('satijalab/seurat-wrappers')
library(reticulate)
library(SCopeLoomR)
library(velocyto.R)

## python3 & scVelo settings
python3 <- Sys.which(names = c("python3.6", "python3"))
python3 <- unname(obj = Filter(f = nchar, x = python3))[1]
reticulate::use_python(python = python3, required = TRUE)
knitr::opts_chunk$set(
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 95),
  message = FALSE,
  warning = FALSE,
  fig.height = 10,
  fig.width = 16
)
## Run in terminal : pip install scvelo --upgrade --quiet
### scVelo original vignette - https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=xc-iI426LGMt

### If there's no loom file, should make it first with velocyto - https://velocyto.org/velocyto.py/tutorial/cli.html#running-velocyto
#options(buildtools.check = function(action) TRUE )
#library(devtools)
#install_github("ssun1116/velocyto.R")
# Not required anymore (maybe) :  openmp & clang4 download : https://thecoatlessprofessor.com/programming/cpp/openmp-in-r-on-os-x/ - solve 'lboost_filesystem' error

library(velocyto.R)
## Loom data contains of spliced / unspliced RNA
loomdat = SeuratWrappers::ReadVelocity(file = "Data/velocyto/Chromium033.loom", engine = "hdf5r")
sc <- as.Seurat(x = loomdat)
sc <- SCTransform(object = sc, assay = "spliced")

ndb_sc <- gsub("Chromium033:", "", colnames(sc))
ndb_sc <- gsub("x", "-1", ndb_sc)
sc <- RenameCells(sc, new.names = ndb_sc) # Ready for merged into sample data
sc <- subset(sc, cells = colnames(helpless),features = rownames(helpless))

helpless = readRDS("Data/SeuratRDS/chromium033_PC10_res0.5.rds")
DimPlot(helpless, reduction = "umap", label = T, label.size = 5)
# Embeddings information - UMAP Plotting
#umap_helpless = cbind("Barcode" = rownames(Embeddings(object = helpless, reduction = "umap")), Embeddings(object = helpless, reduction = "umap"))

sc@reductions <- helpless@reductions
sc@active.ident <- helpless@active.ident
DimPlot(sc, reduction = "umap", label = T, label.size = 5)
sc <- RunVelocity(object = sc, deltaT = 1, kCells = 25, fit.quantile = 0.02) # done!

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sc)))
names(x = ident.colors) <- levels(x = sc)
cell.colors <- ident.colors[Idents(object = sc)]
names(x = cell.colors) <- colnames(x = sc)
show.velocity.on.embedding.cor(emb = Embeddings(object = sc, reduction = "umap"), vel = Tool(object = sc, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

## Loop

samplelist = c("Chromium033", "Chromium034", "Chromium035", "Chromiu040")
for (sample in samplelist){
  ## Loom data contains of spliced / unspliced RNA
  loomdat = SeuratWrappers::ReadVelocity(file = paste("Data/velocyto/", sample, ".loom", sep = ""), engine = "hdf5r")
  sc <- as.Seurat(x = loomdat)
  sc <- SCTransform(object = sc, assay = "spliced")
  
  sobj = readRDS(paste("Data/SeuratRDS/", sample, "_PC10_res0.5.rds", sep = ""))
  DimPlot(sobj, reduction = "umap", label = T, label.size = 5)
  # Embeddings information - UMAP Plotting
  #umap_helpless = cbind("Barcode" = rownames(Embeddings(object = helpless, reduction = "umap")), Embeddings(object = helpless, reduction = "umap"))
  
  ndb_sc <- gsub(paste(sample, ":", sep = ""), "", colnames(sc))
  ndb_sc <- gsub("x", "-1", ndb_sc)
  sc <- RenameCells(sc, new.names = ndb_sc) # Ready for merged into sample data
  sc <- subset(sc, cells = colnames(sobj),features = rownames(sobj))
  
  sc@reductions <- sobj@reductions
  sc@active.ident <- sobj@active.ident
  DimPlot(sc, reduction = "umap", label = T, label.size = 5)
  sc <- RunVelocity(object = sc, deltaT = 1, kCells = 25, fit.quantile = 0.02) # done!
  
  ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sc)))
  names(x = ident.colors) <- levels(x = sc)
  cell.colors <- ident.colors[Idents(object = sc)]
  names(x = cell.colors) <- colnames(x = sc)
  p <- show.velocity.on.embedding.cor(emb = Embeddings(object = sc, reduction = "umap"), vel = Tool(object = sc, 
                                                                                                    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                                      cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                                      do.par = FALSE, cell.border.alpha = 0.1)
  
  save_plot(paste("Figures/RNAvelocity/velocity_", samplename, "_0823.pdf", sep = ""), last_plot(), base_height =10, base_width = 15)
  
}

#######

library(velocyto.R)
## Loom data contains of spliced / unspliced RNA
loomdat = SeuratWrappers::ReadVelocity(file = "Data/velocyto/Chromium034.loom", engine = "hdf5r")
sc <- as.Seurat(x = loomdat)
sc <- SCTransform(object = sc, assay = "spliced")

coghelp = readRDS("Data/SeuratRDS/chromium034_PC10_res0.5.rds")
DimPlot(coghelp, reduction = "umap", label = T, label.size = 5)
# Embeddings information - UMAP Plotting
#umap_coghelp = cbind("Barcode" = rownames(Embeddings(object = coghelp, reduction = "umap")), Embeddings(object = coghelp, reduction = "umap"))

ndb_sc <- gsub("Chromium034:", "", colnames(sc))
ndb_sc <- gsub("x", "-1", ndb_sc)
sc <- RenameCells(sc, new.names = ndb_sc) # Ready for merged into sample data
sc <- subset(sc, cells = colnames(coghelp),features = rownames(coghelp))

sc@reductions <- coghelp@reductions
sc@active.ident <- coghelp@active.ident
DimPlot(sc, reduction = "umap", label = T, label.size = 5)
sc <- RunVelocity(object = sc, deltaT = 1, kCells = 25, fit.quantile = 0.02) # done!

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sc)))
names(x = ident.colors) <- levels(x = sc)
cell.colors <- ident.colors[Idents(object = sc)]
names(x = cell.colors) <- colnames(x = sc)
show.velocity.on.embedding.cor(emb = Embeddings(object = sc, reduction = "umap"), vel = Tool(object = sc, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

######

ldat <- ReadVelocity(file = "Data/velocyto/Chromium033.loom", engine = "hdf5r")
bm <- as.Seurat(x = ldat)
bm <- SCTransform(object = bm, assay = "spliced")

bm_sc <- gsub("Chromium033:", "", colnames(bm))
bm_sc <- gsub("x", "-1", bm_sc)
bm <- RenameCells(bm, new.names = bm_sc) # Ready for merged into sample data
bm <- subset(bm, cells = colnames(helpless),features = rownames(helpless))

bm <- RunPCA(object = bm, verbose = FALSE)
bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

