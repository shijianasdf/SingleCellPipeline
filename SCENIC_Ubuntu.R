options(stringsAsFactors = F)
library(tidyverse)
library(dplyr)
library(SCENIC) 
library(AUCell)
library(RcisTarget)
library(SCopeLoomR) # install.packages("hdf5r", configure.args = c("--with-hdf5=/opt/homebrew/Cellar/hdf5/1.12.1/bin/h5cc"))
library(KernSmooth)
library(BiocParallel)
library(ggplot2)
library(data.table) ## may need re download.
library(grid)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(Seurat)
options(timeout = max(3000, getOption("timeout"))) ## download.file timeout

### Species-specific database
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs

for(featherURL in dbFiles){
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}

##################### HARMONY #########################
## Get data from sce object:
sce = readRDS("Data/UMAP.Harmony_cluster_0to16_per100_with_naive_0714.rds")
exprMat <- sce@assays$RNA@counts %>% as.matrix()
cellInfo <- data.frame(seuratCluster=Idents(sce))

## Saving into loom
setwd("/data2/scRNA_seq_SYK")
dir.create("Data/SCENIC_loom.data")
loom <- build_loom("Data/SCENIC_loom.data/Harmony_0813.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

### Load data
loomPath = "Data/SCENIC_loom.data/Harmony_0813.loom"
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)
#saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo"))

### Initialize settings
scenicOptions <- initializeScenic(org="mgi", dbDir="Resources", nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# # Color to assign to the variables (same format as for NMF::aheatmap)
# colVars <- list(CellType=c("microglia"="forestgreen", 
#                            "endothelial-mural"="darkorange", 
#                            "astrocytes_ependymal"="magenta4", 
#                            "oligodendrocytes"="hotpink", 
#                            "interneurons"="red3", 
#                            "pyramidal CA1"="skyblue", 
#                            "pyramidal SS"="darkblue"))
# colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
# saveRDS(colVars, file="int/colVars.Rds")
# #plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))
# # only legend..!!

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)


### Build and score the GRN
#exprMat_log <- log2(exprMat+1) # why not using the filtered one...??
scenicOptions <- readRDS("int/scenicOptions.Rds")
#scenicOptions@settings$verbose <- TRUE
#scenicOptions@settings$nCores <- 10
#scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
#scenicOptions@settings$nCores <- 1
#(https://github.com/aertslab/AUCell/blob/master/R/02_calcAUC.R) # use this function..
#source('02_calcAUC.R')
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
#scenicOptions@settings$nCores <- 10
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Optional: Binarize activity
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered_log)
savedSelections <- shiny::runApp(aucellApp)
# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions, exprMat = exprMat_filtered_log)
#tsneAUC(scenicOptions, aucType="AUC") # choose settings
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


# Clustering / dimensionality reduction on the regulon activity
#nPcs <- c(5) # For toy dataset
nPcs <- c(5,15,50)
#scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

# to view
pdf(file = 'output/CellType_1.pdf', width = 8, height = 8)
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=F, varName='CellType', cex=.5)
dev.off()

# Using only "high-confidence" regulons (normally similar)
pdf(file = 'output/CellType_2.pdf', width = 8, height = 8)
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)
dev.off()

# Using fixed combinations
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
library(SCopeLoomR)
scenicOptions@fileNames$output["loomFile",] <- "output/mouseBrain_SCENIC.loom"
export2loom(scenicOptions, exprMat)
# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org
# Projection the AUC and TF expression onto t-SNEs
#exprMat_log <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered_log) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)
print(tsneFileName(scenicOptions))

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
pdf(file = 'output/TF_expression.pdf', width = 10, height = 8)
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_filtered_log, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Dlx5", "Sox10", "Sox9","Irf1", "Stat6")],], plots="Expression")
dev.off()

# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_filtered_log, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

# density plot
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
Cairo::CairoPDF("output/density_plot.pdf", width=6, height=6)
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)
dev.off()

#Show several regulons simultaneously:
#par(bg = "black")
Cairo::CairoPDF("output/regulons_plot.pdf", width=8, height=4)
par(mfrow=c(1,2))
regulonNames <- c( "Dlx5","Sox10")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
regulonNames <- list(red=c("Sox10", "Sox8"),
                     green=c("Irf1"),
                     blue=c( "Tef"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
dev.off()

## GRN: Regulon targets and motifs
# Genes included in the regulons:
regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Dlx5", "Irf1")]
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

# output/Step2_MotifEnrichment_preview.html in detail/subset:
# getOutName(scenicOptions, "s2_motifEnrichment")
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=10))


## Regulators for known cell types or clusters
# Average Regulon Activity by cluster
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
#
pdf(file = 'output/regulonActivity_byCellType_heatmap.pdf', width = 5, height = 4)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
dev.off()

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

## Binarized version
minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pdf(file = 'output/regulonActivity_byCellType_binarized_heatmap.pdf', width = 5, height = 4)
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))
dev.off()

topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)

# Cell-type specific regulators (RSS): (Useful for big analysis with many cell types, to identify the cell-type specific regulons.)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
