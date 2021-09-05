# SNU_scRNA_Immune
Organize the workflow of myeloid cell single cell rna-seq analysis. 
## Running cellRanger count
- Cell Ranger is a set of analysis pipelines that process Chromium single-cell data to align reads, generate feature-barcode matrices, perform clustering and other secondary analysis, and more.
- `cellranger count` takes FASTQ files from cellranger mkfastq and performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis. The count pipeline can take input from multiple sequencing runs on the same GEM well. cellranger count also processes Feature Barcode data alongside Gene Expression reads.
3. Seurat
4. Harmony
5. Monocle - Pseudotime trajectory analysis
+ FGSEA. Fisher Test. DEG ...
