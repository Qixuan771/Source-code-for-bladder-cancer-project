mkdir raw_bw_pool
computeMatrix scale-regions -p 22 \
-S $TF_bw/RT4_FOXA1_pool*.bw $TF_bw/RT4_GATA3*.bw $ATAC_bw/RT4*.bw $K27ac_bw/RT4*.bw \
-R *enh*.txt \
--beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
--skipZeros \
-o raw_bw_pool/RT4_ATAC_enh_wK27ac_wTF_enh.gz

plotHeatmap -m raw_bw_pool/RT4_ATAC_enh_wK27ac_wTF_enh.gz \
-out raw_bw_pool/RT4_ATAC_enh_wK27ac_wTF_enh.pdf \
--dpi 300 \
--samplesLabel RT4_FOXA1 RT4_GATA3 RT4_ATAC RT4_H3K27ac \
--whatToShow 'heatmap and colorbar' \
--zMin 0 --zMax 40 \
--colorList 'white,purple' 'white,magenta' 'white,cyan' 'white,orange' 
