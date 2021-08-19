
## Code for installing packages needed for SCENIC analysis

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## Required
BiocManager::install(c("AUCell", "RcisTarget"), force = TRUE)
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"), force = TRUE)
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"), force = TRUE)
# To support parallel execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"), force = TRUE)

## Couldn't download these ... Try again when needed
# # To export/visualize in http://scope.aertslab.org
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("hhoeflin/hdf5r", configure.args = c("--with-hdf5=/opt/homebrew/bin/h5cc"))
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC") # 1.2.4

#install.packages(c("jsonlite", "evaluate"))