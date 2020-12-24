f=${HT_bam[$i]};
    
samtools index $f
    
total=$(samtools view -c $f)
    
reads_in_peaks=$(bedtools sort -i ${HT[$i]}|\
                 bedtools merge -i stdin|\
                 bedtools intersect -u -nonamecheck -a $f -b stdin -ubam|\
                 samtools view -c)
    
calc(){ awk "BEGIN { print "$*" }"; }

FRiP=$(calc ${reads_in_peaks}/${total})
    
echo ___For HT1376 - ${f} ___
echo total: $total
echo reads_in_peaks: $reads_in_peaks
echo FRiP: $FRiP
    
scalefactor=$(calc 10000000/${reads_in_peaks})
    
bamCoverage \
-b $f \
-o RiPNorm/${f%.*}_RiPnorm.bw \
--scaleFactor $scalefactor \
--numberOfProcessors 22;
