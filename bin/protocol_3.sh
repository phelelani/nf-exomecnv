#!/bin/bash

VCF=$1

## 1. FILTER CNV BASED ON FREQUENCY ESTIMATED FROM THE SAMPLE
##-- SCRIPT TO CONVERT XHMM CNV CALLS TO PLINK FORMAT (REQUIRES PLINK 1.7)
grep "^#CHROM" $VCF | awk '{for (i = 10; i <= NF; i++) print $i,1,0,0,1,1}' > DATA.fam
xcnv_to_cnv DATA.xcnv > DATA.cnv
/home/phelelani/applications/plink/plink --noweb --cfile DATA --cnv-make-map --out DATA

##-- FIND ANY CNV THAT OVERLAPS 50% OF ANOTHER CNV THAT OCCURS IN MORE THAN 10% OF ALL SAMPLES
##-- (A RELATIVELY LIBERAL FREQUENCY THRESHOLD TO REMOVE ONLY COMMON CNV AND ARTIFACTS).
## CALCULATE THE NUMBER OF SAMPLES IN THE DATA
NUM_SAMPLES=`cat DATA.fam | wc -l`

## SET THE FREQUENCY THRESHOLD AT 10% OF THE TOTAL NUMBER OF SAMPLES
THRESH_NUM=`echo "0.1 * $NUM_SAMPLES" | bc | awk '{x = $1; if (x != int(x)) {x = int(x)+1} print x}'`
/home/phelelani/applications/plink/plink --noweb --cfile DATA --cnv-freq-exclude-above $THRESH_NUM --cnv-overlap 0.5 --cnv-write --out DATA.maf_0.1
cat DATA.maf_0.1.cnv | awk '{if (NR>1) print $3":"$4"-"$5}' | sort | uniq | sed 's/^/chr/g' > DATA.maf_0.1.CNV_regions.txt

## EXCLUDE ALL SUCH COMMON REGIONS FROM THE VCF FILE
cat $VCF | awk '{if (substr($1,1,1)=="#") {print $_} else {found=-1; cmd = "grep -w "$3" DATA.maf_0.1.CNV_regions.txt"; cmd | getline found; close(cmd);
if (found != -1) {print $_}}}' > DATA.maf_0.1.vcf

## 2. DETECT PUTATIVE DE NOVO CNV BY STRICT USE OF XHMM QUALITY SCORES
## CREATE A NEW PLINK/SEQ PROJECT THAT WILL CONTAIN THE XHMM-CALLED VARIANTS:
/home/phelelani/applications/plinkseq/pseq DATA new-project

## LOAD A PLINK/SEQ PEDINFO FILE CONTAINING THE CHILD-PARENT PEDIGREE RELATIONSHIPS (NAMED ‘DATA.PEDINFO’ HERE):
/home/phelelani/applications/plinkseq/pseq DATA load-pedigree --file DATA.pedinfo

## LOAD THE FREQUENCY-FILTERED VCF FILE INTO THE PLINK/SEQ PROJECT:
/home/phelelani/applications/plinkseq/pseq DATA load-vcf --vcf DATA.maf_0.1.vcf

## NOW, ASSUMING YOU'VE CHOSEN A QUALITY SCORE FILTER OF 60 TO CALL DE NOVO CNV ON THE AUTOSOMES, RUN:
/home/phelelani/applications/plinkseq/pseq DATA cnv-denovo --mask reg.ex=chrX,chrY --minSQ 60 --minNQ 60 --out DATA



## PROTOCOL 2
pseq . loc-intersect --group refseq --locdb locdb --file EXOME.interval_list --out annotated_targets.refseq
