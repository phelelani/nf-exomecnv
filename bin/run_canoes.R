#!/usr/bin/env Rscript

## LOAD THE REQUIRED PACKAGES - INSTALL IF MISSING!
if(!require(nnls)){
    install.packages("nnls", dependencies=TRUE)
    library("nnls")
}
if(!require(Hmisc)){
    install.packages("Hmisc", dependencies=TRUE)
    library("Hmisc")
}
if(!require(mgcv)){
    install.packages("mgcv", dependencies=TRUE)
    library("mgcv")
}
if(!require(plyr)){
    install.packages("plyr", dependencies=TRUE)
    library("plyr")
}

## GET USER INPUT
args = commandArgs(trailingOnly=TRUE)
in.gc <- args[1]
in.reads <- args[2]
in.samples <- args[3]

## READ IN DATA
chr <- gsub("_canoes_reads_new.txt", "", in.reads)
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
xcnv.pass <- vector('list', length(sample.names))
xcnv.fail <- vector('list', length(sample.names))

## CALL CNVS IN EACH SAMPLE
for (i in 1:length(sample.names)){
    cat("\nI am now finding CNVs for sample: ", sample.names[i], "\n")
    tryCatch(
    {
        xcnv.pass[[i]] <- CallCNVs(sample.names[i], canoes.reads)
        if ( nrow(xcnv.pass[[i]]) >= 50 ) {
            cat("-- This sample ", sample.names[i], " has more than 50 CNV's - I'm removing it!\n")
            xcnv.fail[[i]] <- xcnv.pass[[i]]
            xcnv.pass[[i]] <- NULL
        }
    },
    error = function(e) {
        message('An Error Occurred')
        print(e)
    },
    warning = function(w) {
        message('A Warning Occurred')
        print(w)
    }
    )
    cat("-- Done finding CNVs\n")
}

## COMBINE THE RESULTS INTO ONE DATA FRAME - WRITE TO CSV AND SAVE THE SESSION
xcnvs.pass <- do.call('rbind', xcnv.pass)
xcnvs.fail <- do.call('rbind', xcnv.fail)
write.table(xcnvs.pass, file=paste0(chr, "_CNVs_pass.csv"), sep = "\t", row.names=FALSE, col.names = TRUE, quote=FALSE)
write.table(xcnvs.fail, file=paste0(chr, "_CNVs_fail.csv"), sep = "\t", row.names=FALSE, col.names = TRUE, quote=FALSE)
save.image(file=paste0(chr, "_xcnvs.RData"))

## CREATE DIRECTORIES TO PUT GENOTYPE INFO AND PLOTS
dir_plots <- paste0(chr, "_CNV_plots")
dir_geno <- paste0(chr, "_CNV_genotype")
dir.create(dir_plots)
dir.create(dir_geno)

for (the.sample in unique(xcnvs.pass$SAMPLE)){
    cat("\nI am now genotyping and plotting CNVs for the sample: ", the.sample,"\n")
    ## GENOTYPE ALL THE CNV CALLS MADE ABOVE PER SAMPLE
    the.name <- paste0(the.sample, "_", chr, "_genotype")
    assign(the.name, GenotypeCNVs(xcnvs.pass, the.sample, canoes.reads))
    write.table(eval(parse(text=the.name)), file=paste0(dir_geno, "/", the.name, ".csv"), sep = "\t", row.names=FALSE,col.names = TRUE,quote=FALSE)
    cat("-- Done genotyping\n")
    ## PLOT ALL THE CNV CALLS PER SAMPLE    
    pdf(paste0(dir_plots, "/", the.sample, "_", chr, "_CNVs_plot.pdf"))
    the.sample.cnvs <- xcnvs.pass[xcnvs.pass$SAMPLE==the.sample,]
    if (the.sample %in% the.sample.cnvs$SAMPLE){
        for (i in 1:nrow(the.sample.cnvs)){
            PlotCNV(canoes.reads, the.sample.cnvs[i, "SAMPLE"], the.sample.cnvs[i, "TARGETS"])        
        }
    }
    dev.off()
    cat("-- Done plotting\n")
}
