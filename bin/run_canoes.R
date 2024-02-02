#!/usr/bin/env Rscript

library("nnls")
library("Hmisc")
library("mgcv")
library("plyr")

## GET USER INPUT
args = commandArgs(trailingOnly=TRUE)
in.gc <- args[1]
in.reads <- args[2]
in.samples <- args[3]

## READ IN DATA
gc <- read.table(in.gc)$V2
canoes.reads <- read.table(in.reads)

## RENAME CANOES.READS COLUMNS
sample.names <- scan(in.samples, what="", sep="\n")
names(canoes.reads) <- c("chromosome", "start", "end", sample.names)

## CREATE A VECTOR OF CONSECUTIVE TARGET IDS
target <- seq(1, nrow(canoes.reads))

## COMBINE THE DATA INTO ONE DATA FRAME
canoes.reads <- cbind(target, gc, canoes.reads)

## CALL CNVS IN THE SAMPLES
source("/home/phelelani/nf-workflows/nf-exomecnv/bin/CANOES.R")

## CREATE A VECTOR TO HOLD THE RESULTS FOR EACH SAMPLE
xcnv.list <- vector('list', length(sample.names))

## CALL CNVS IN EACH SAMPLE
for (i in 1:length(sample.names)){
    xcnv.list[[i]] <- CallCNVs(sample.names[i], canoes.reads)
}

## COMBINE THE RESULTS INTO ONE DATA FRAME
xcnvs <- do.call('rbind', xcnv.list)
write.table(xcnvs, file="Sample_CNVs.csv", sep = "\t", row.names=FALSE, col.names = TRUE, quote=FALSE)

dir.create("Sample_CNV_plots")
dir.create("Sample_Genotyping")
           
for (the.sample in sample.names){
    ## GENOTYPE ALL THE CNV CALLS MADE ABOVE PER SAMPLE
    the.name <- paste0(the.sample, "_genotyping")
    assign(the.name, GenotypeCNVs(xcnvs, the.sample, canoes.reads))
    write.table(eval(parse(text=the.name)), file=paste0("Sample_Genotyping/", the.name, ".csv"), sep = "\t", row.names=FALSE, col.names = TRUE, quote=FALSE)
    
    ## PLOT ALL THE CNV CALLS PER SAMPLE    
    pdf(paste0("Sample_CNV_plots/", the.sample, "_CNVplot.pdf"))
    tmp.xcnvs <- xcnvs[xcnvs$SAMPLE==the.sample,]

    if (the.sample %in% tmp.xcnvs$SAMPLE){
        for (i in 1:nrow(tmp.xcnvs)){
            PlotCNV(canoes.reads, tmp.xcnvs[i, "SAMPLE"], tmp.xcnvs[i, "TARGETS"])        
        }
    }
    dev.off()
}
