#!/bin/bash

NUM_SAMPLES=`ls *.norm.cov.bed | wc -l | awk '{print $1}'`
NUM_WINDOWS=`ls *.norm.cov.bed | head -n 1 | xargs awk '$1 != "X" && $1 != "Y" && $NF == 0 {x++} END {print x}'`

echo -e "$NUM_SAMPLES\t$NUM_WINDOWS" > matrix.txt

ls *.norm.cov.bed | while read FILE
do
    awk '$1 != "X" && $1 != "Y" && $NF != 0 { print $4 }' $FILE | gawk -f $CLAMMS_DIR/transpose.gawk >> matrix.txt
done

svd -d 4 -o svd-output -r dt matrix.txt
ls *.norm.cov.bed | cut -d '.' -f 1 > sample.names.txt
tail -n +2 svd-output-Ut | tr ' ' '\t' | gawk -f $CLAMMS_DIR/transpose.gawk | paste sample.names.txt - > pca.coordinates.txt
