# SNU_scRNA_Immune

### Description
Organize the workflow of myeloid cell single cell rna-seq analysis pipeline. 

## Running `cellRanger count`
- Cell Ranger is a set of analysis pipelines that process Chromium single-cell data to align reads, generate feature-barcode matrices, perform clustering and other secondary analysis, and more.
- `cellranger count` takes FASTQ files from cellranger mkfastq and performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis. The count pipeline can take input from multiple sequencing runs on the same GEM well. cellranger count also processes Feature Barcode data alongside Gene Expression reads.
## Analysis with **Seurat**
- main page : https://satijalab.org/seurat/index.html
- Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat enables to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data.
4. Harmony
5. Monocle - Pseudotime trajectory analysis
+ FGSEA. Fisher Test. DEG ...
