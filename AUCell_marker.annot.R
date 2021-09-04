library(Seurat)
library(msigdbi)
library(tidyverse)
library(GSEABase)
library(AUCell)

harmony = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")

######   Marker-based automatic annotation  ######

## SCINA
# Import the marker genes as a GMT file and store as a variable
markers <- msigdbi::read.gmt('Resources/c8.all.v7.2.symbols.gmt')
# Convert the expression data from Seurat object into a matrix data structure
exprMatrix <- data.frame(Seurat::GetAssayData(harmony))
exprMatrix2 <- as.matrix(Seurat::GetAssayData(harmony))
gg <- rownames(exprMatrix)
gg <- gg %>% toupper
rownames(exprMatrix) <- gg
col <- colnames(exprMatrix2)
colnames(exprMatrix) <- col
exprMatrix = as.matrix(exprMatrix)
# Run SCINA on the query data using the marker genes to identify cell types
# Specifying rm_overlap = FALSE allows the same marker gene to specify multiple cell types which
# may be useful if identifying cell subtypes or other similar types of cells.
# Specifying allow_unknown = TRUE allows cells to be labeled as "unknown" instead of being assigned a low-confident label.
predictions.scina = SCINA::SCINA(exp = exprMatrix, signatures = markers$genesets,
                                 rm_overlap = FALSE, allow_unknown = TRUE) ## Not applicable for large geneset, cannot process.
# Add SCINA annotation information to each cell in Seurat object
colData(query_sce)$SCINA <- predictions.scina$cell_labels

# Make a UMAP and add the SCINA cell-type annotations
scater::plotReducedDim(query_sce, dimred="UMAP", colour_by="SCINA") +
  ggplot2::theme(legend.position = "bottom",
                 legend.text = ggplot2::element_text(size = 4))


## AUCell
# Setup
# To support paralell execution:
BiocManager::install(c("doMC", "doRNG","doSNOW"), force = TRUE)
# For the main example:
BiocManager::install(c("mixtools", "GEOquery", "SummarizedExperiment"), force = TRUE)
# For the examples in the follow-up section of the tutorial:
BiocManager::install(c("DT", "plotly", "NMF", "shiny", "rbokeh", "dynamicTreeCut","R2HTML","Rtsne", "zoo"), force = TRUE)
devtools::install_github("talgalili/d3heatmap")
library(AUCell)
dir.create("AUCell")
setwd("AUCell") 

# Convert the expression data from Seurat object into a matrix data structure
exprMatrix <- data.frame(Seurat::GetAssayData(harmony))
exprMatrix2 <- as.matrix(Seurat::GetAssayData(harmony))
gg <- rownames(exprMatrix)
gg <- gg %>% toupper
rownames(exprMatrix) <- gg
col <- colnames(exprMatrix2)
colnames(exprMatrix) <- col
exprMatrix = as.matrix(exprMatrix) # Modify rowname and colnames for appropriate format.
# Load Geneset
geneSets <- getGmt('../Resources/c8.all.v7.2.symbols.gmt')
## Filter Bone Marrow geneset
names <- names(geneSets)
names = names[grep('BONE_MARROW', names)]
geneSets = geneSets[names]
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))

# 1. Build gene-expression rankings for each cell
cells_rankings <- AUCell_buildRankings(exprMatrix)
#save(cells_rankings, file="cells_rankings.RData")
#cells_rankings = load("cells_rankings.RData")

# 2. Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="cells_AUC.RData")
#cells_AUC = load("cells_AUC.RData")

# 3. Determine the cells with the given gene signatures or active gene sets
set.seed(123)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

DimPlot(harmony, reduction = "umap", label = TRUE, pt.size = .1)
har.umap = harmony@reductions$umap@cell.embeddings
plot(har.umap, pch = 16, cex = .3)

selectedThresholds <- getThresholdSelected(cells_assignment)
for(geneSetName in names(selectedThresholds)){
  pdf(file = paste("UMAP.plot.", geneSetName, "_AUC.scored.Harmony.0831.pdf", sep = ""), width = 8, height = 8)
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0){
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    # Plot
    plot(har.umap, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(har.umap)], pch=16) 
  }
  dev.off()
}


## AUCell analysis 2 - cellmarker geneset

rm(list = ls())
setwd("AUCell")
harmony = readRDS("../Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")

# Convert the expression data from Seurat object into a matrix data structure
exprMatrix <- as.matrix(Seurat::GetAssayData(harmony))
# exprMatrix2 <- as.matrix(Seurat::GetAssayData(harmony))
# gg <- rownames(exprMatrix)
# gg <- gg %>% toupper
# rownames(exprMatrix) <- gg
# col <- colnames(exprMatrix2)
# colnames(exprMatrix) <- col
# exprMatrix = as.matrix(exprMatrix) # Modify rowname and colnames for appropriate format.

# Load Geneset
geneset1 = read.delim("../Resources/References/Mouse_cell_markers.txt")
geneset2 = read.delim("../Resources/References/Single_cell_markers.txt")
geneset = rbind.data.frame(geneset1, geneset2) %>% filter(speciesType == "Mouse" & tissueType %in% c("Bone marrow", "Spleen", "Lymph node") &
                                                            cancerType == "Normal") %>% unique
geneset_ls = list()
for (cellname in unique(geneset$cellName)){
  d = geneset %>% filter(cellName  == cellname)
  str = character()
  for (i in 1:nrow(d))
    str = paste(str, d[i,]$geneSymbol, by = ",")
  str = gsub(" ", "", str)
  str = unlist(strsplit(str, ",")) %>% unique
  geneset_ls[[cellname]] = GeneSet(geneIds=str, setName = cellname)
}
geneSets <- GeneSetCollection(geneset_ls)
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
save(geneSets, file = "geneSets_cellmarker.RData")
load("geneSets_cellmarker.RData")

# 1. Build gene-expression rankings for each cell
cells_rankings <- AUCell_buildRankings(exprMatrix)
save(cells_rankings, file="cells_rankings.RData")
#cells_rankings = load("cells_rankings.RData")

# 2. Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="cells_AUC.cellmarker.RData")
#cells_AUC = load("cells_AUC.RData")

# 3. Determine the cells with the given gene signatures or active gene sets
set.seed(123)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

DimPlot(harmony, reduction = "umap", label = TRUE, pt.size = .1)
har.umap = harmony@reductions$umap@cell.embeddings
plot(har.umap, pch = 16, cex = .3)

selectedThresholds <- getThresholdSelected(cells_assignment)
for(geneSetName in names(selectedThresholds)){
  pdf(file = paste("UMAP.plot.cellmarker.", geneSetName, "_AUC.scored.Harmony.0831.pdf", sep = ""), width = 8, height = 8)
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0){
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    # Plot
    plot(har.umap, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(har.umap)], pch=16) 
  }
  dev.off()
}


## AUCell analysis 3 - DC reference gene list

rm(list = ls())
dir.create("AUCell")
setwd("AUCell")
harmony = readRDS("../Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")

# Convert the expression data from Seurat object into a matrix data structure
exprMatrix <- as.matrix(Seurat::GetAssayData(harmony))
# exprMatrix2 <- as.matrix(Seurat::GetAssayData(harmony))
# gg <- rownames(exprMatrix)
# gg <- gg %>% toupper
# rownames(exprMatrix) <- gg
# col <- colnames(exprMatrix2)
# colnames(exprMatrix) <- col
# exprMatrix = as.matrix(exprMatrix) # Modify rowname and colnames for appropriate format.

# Load Geneset
## Phenotype encoding genes
cdc1_marker <- c("Cd8a", "Xcr1", "Cd24", "Clec9a", "Ly75", "Zbtb46", "Itgae", "Irf8", "Spi1", "Batf3", "Id2")
cdc1_marker <- str_split(cdc1_marker, ", ") %>% unlist # char type.
cdc2_marker <- c("Itgam", "Sirpa", "Cd4", "Clec10a", "Clec12a", "Zbtb46", "Irf4", "Spi1")
cdc2_marker <- str_split(cdc2_marker, ", ") %>% unlist # char type.
mo_marker <- c("Itgam", "Cd209", "Mrc1", "Mrc2", "Fcgr1", "Adgre1", "Sirpa", "Cd14", "Irf4", "Spi1")
mo_marker <- str_split(mo_marker, ", ") %>% unlist # char type.
pdc_marker <- c("Ptprc", "Bst2", "Ly49q", "Siglech", "Ccr9", "Tcf4", "Zeb2")
pdc_marker <- str_split(pdc_marker, ", ") %>% unlist # char type.

marker_ls = list()
marker_ls[["cDC1"]] = GeneSet(geneIds=cdc1_marker, setName = "cDC1")
marker_ls[["cDC2"]] = GeneSet(geneIds=cdc2_marker, setName = "cDC2")
marker_ls[["moDC"]] = GeneSet(geneIds=mo_marker, setName = "moDC")
marker_ls[["pDC"]] = GeneSet(geneIds=pdc_marker, setName = "pDC")

geneSets <- GeneSetCollection(marker_ls)
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
save(geneSets, file = "geneSets_DCmarker.RData")
#load("geneSets_cellmarker.RData")

# 1. Build gene-expression rankings for each cell
cells_rankings <- AUCell_buildRankings(exprMatrix)
save(cells_rankings, file="cells_rankings.RData")
#cells_rankings = load("cells_rankings.RData")

# 2. Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="cells_AUC.cellmarker.RData")
#cells_AUC = load("cells_AUC.RData")

# 3. Determine the cells with the given gene signatures or active gene sets
set.seed(123)
par(mfrow=c(2,2)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

DimPlot(harmony, reduction = "umap", label = TRUE, pt.size = .1)
har.umap = harmony@reductions$umap@cell.embeddings
plot(har.umap, pch = 16, cex = .3)

selectedThresholds <- getThresholdSelected(cells_assignment)
for(geneSetName in names(selectedThresholds)){
  pdf(file = paste("UMAP.plot.markerlist.", geneSetName, "_AUC.scored.Harmony.0902.pdf", sep = ""), width = 8, height = 8)
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0){
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    # Plot
    plot(har.umap, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(har.umap)], pch=16) 
  }
  dev.off()
}



## AUCell analysis 3 - DC reference gene list

rm(list = ls())
dir.create("AUCell")
setwd("AUCell")
harmony = readRDS("../Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")

# Convert the expression data from Seurat object into a matrix data structure
exprMatrix <- as.matrix(Seurat::GetAssayData(harmony))
# exprMatrix2 <- as.matrix(Seurat::GetAssayData(harmony))
# gg <- rownames(exprMatrix)
# gg <- gg %>% toupper
# rownames(exprMatrix) <- gg
# col <- colnames(exprMatrix2)
# colnames(exprMatrix) <- col
# exprMatrix = as.matrix(exprMatrix) # Modify rowname and colnames for appropriate format.

# Load Geneset
## Phenotype encoding genes
cdc1_marker <- c("Cd8a", "Xcr1", "Cd24", "Clec9a", "Ly75", "Zbtb46", "Itgae", "Irf8", "Spi1", "Batf3", "Id2")
cdc1_marker <- str_split(cdc1_marker, ", ") %>% unlist # char type.
cdc2_marker <- c("Itgam", "Sirpa", "Cd4", "Clec10a", "Clec12a", "Zbtb46", "Irf4", "Spi1")
cdc2_marker <- str_split(cdc2_marker, ", ") %>% unlist # char type.
mo_marker <- c("Itgam", "Cd209", "Mrc1", "Mrc2", "Fcgr1", "Adgre1", "Sirpa", "Cd14", "Irf4", "Spi1")
mo_marker <- str_split(mo_marker, ", ") %>% unlist # char type.
pdc_marker <- c("Ptprc", "Bst2", "Ly49q", "Siglech", "Ccr9", "Tcf4", "Zeb2")
pdc_marker <- str_split(pdc_marker, ", ") %>% unlist # char type.

marker_ls = list()
marker_ls[["cDC1"]] = GeneSet(geneIds=cdc1_marker, setName = "cDC1")
marker_ls[["cDC2"]] = GeneSet(geneIds=cdc2_marker, setName = "cDC2")
marker_ls[["moDC"]] = GeneSet(geneIds=mo_marker, setName = "moDC")
marker_ls[["pDC"]] = GeneSet(geneIds=pdc_marker, setName = "pDC")

geneSets <- GeneSetCollection(marker_ls)
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
save(geneSets, file = "geneSets_DCmarker.RData")
#load("geneSets_cellmarker.RData")

# 1. Build gene-expression rankings for each cell
cells_rankings <- AUCell_buildRankings(exprMatrix)
save(cells_rankings, file="cells_rankings.RData")
#cells_rankings = load("cells_rankings.RData")

# 2. Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="cells_AUC.cellmarker.RData")
#cells_AUC = load("cells_AUC.RData")

# 3. Determine the cells with the given gene signatures or active gene sets
set.seed(123)
par(mfrow=c(2,2)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

DimPlot(harmony, reduction = "umap", label = TRUE, pt.size = .1)
har.umap = harmony@reductions$umap@cell.embeddings
plot(har.umap, pch = 16, cex = .3)

selectedThresholds <- getThresholdSelected(cells_assignment)
for(geneSetName in names(selectedThresholds)){
  pdf(file = paste("UMAP.plot.markerlist.", geneSetName, "_AUC.scored.Harmony.0902.pdf", sep = ""), width = 8, height = 8)
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0){
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    # Plot
    plot(har.umap, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(har.umap)], pch=16) 
  }
  dev.off()
}








