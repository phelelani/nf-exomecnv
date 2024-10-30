#!/bin/bash

## GET FILE HEADERS
grep --no-filename --no-group-separator -A1 "## METRICS CLASS" *.hs_metrics.txt | grep -v "## METRICS CLASS" | uniq | awk -v OFS="\t" '{ print "SAMPLE",$9,$32,$48,$52,$63,$64 }' > hs_metrics.txt
grep --no-filename --no-group-separator -A1 "## METRICS CLASS" *.insert_size_metrics.txt | grep -v "## METRICS CLASS" | uniq | awk  -v OFS="\t" '{ print "SAMPLE",$6}' > insert_size_metrics.txt

## PRINT OUT FOR HS METRICS
for hs_metrics in *.hs_metrics.txt
do
    sample="${hs_metrics%.hs_metrics.txt}"
    sed -n '/^BAIT_SET\(.*\)$/,/^$/p' $hs_metrics |\
        sed '1d;$d' | awk -v OFS="\t" -v sample="$sample" '{ print sample,$9,$32,$48,$52,$63,$64 }'
done >> hs_metrics.txt

## PRINT OUT FOR INSERT SIZE METRICS
for ins_metrics in *.insert_size_metrics.txt
do
    sample="${ins_metrics%.insert_size_metrics.txt}"
    sed -n '/^MEDIAN_INSERT_SIZE\(.*\)$/,/^$/p' $ins_metrics |\
        sed '1d;$d' | awk -v OFS="\t" -v sample="$sample" '{ print sample,$6}'
done >> insert_size_metrics.txt

join -j 1 -t '	' hs_metrics.txt insert_size_metrics.txt > qcs_metrics
