## Analyze RNA velocity quantifications stored in a Seurat object using scVelo.
options(stringsAsFactors = F)
library(dplyr)
library(tidyverse)
library(Seurat)
library(SeuratDisk) # remotes::install_github("mojaveazure/seurat-disk")
library(SeuratWrappers) # remotes::install_github('satijalab/seurat-wrappers')
library(reticulate)

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
### scVelo original vignette link
### : https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=xc-iI426LGMt


### If there's no loom file, should make it first with velocyto 
### : https://velocyto.org/velocyto.py/tutorial/cli.html#running-velocyto
#BiocManager::install("pcaMethods", force = TRUE)
library(devtools)
install_github("velocyto-team/velocyto.R")
loomdat = ReadVelocity(file = "Data/velocyto/Chromium033.loom")
seuratobject = readRDS("Data/SeuratRDS/chromium033_PC10_res0.5.rds")
