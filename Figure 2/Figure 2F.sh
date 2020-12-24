
inp=$1
out=200107_ATACSeq_RiPNorm_TCGA_multiBigwigSummary_pool_${inp}

multiBigwigSummary BED-file \
--BED ${inp} \
-b ${ATAC}/*.bw ${TCGA}/*.bw \
-o ${out}.npz \
--outRawCounts ${out}.txt

plotCorrelation \
-in ${out}.npz \
--corMethod pearson --skipZeros --removeOutliers \
--plotTitle "Pearson Correlation of BLCA ATACSeqs" \
--whatToPlot heatmap --colorMap bwr --plotNumbers \
-o ${out}_PearsonCorr_hmap_bigwigScores.pdf   \
--outFileCorMatrix ${out}_PearsonCorr_hmap_bigwigScores.tab

