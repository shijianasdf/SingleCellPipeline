declare -a arr=("Chromium033" "Chromium034" "Chromium035" "Chromium040");
for sample in "${arr[@]}"; do echo "$i" ; 
#mkdir ${sample}/velocyto_umap ; 
(velocyto run \
-b ${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o ${sample}/velocyto_umap \
${sample}/outs/possorted_genome_bam.bam /data2/scRNA_seq/Resources/refdata-gex-mm10-2020-A/genes/genes.gtf &) ; 
done
