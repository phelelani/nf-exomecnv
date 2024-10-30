#!/usr/bin/env Rscript

library(ggplot2)
library(FNN)

## GET USER INPUT
args = commandArgs(trailingOnly=TRUE)
in.pca <- args[1]
in.qcs <- args[2]
in.sex <- args[3]

## GET THE PCA DATA
pca.data <- read.table(in.pca, col.names=c("SAMPLE", "PC1", "PC2", "PC3", "PC4"), colClasses=c("character", rep("numeric", 4)))

## PLOT PCA
pdf("pca.pdf")
ggplot(pca.data, aes(x = PC1, y = PC2)) +
    geom_point()
dev.off()

## load("/home/phelelani/applications/clamms/data/example_qcs.Rdata")
qcs.data <- read.table(in.qcs, header=TRUE, sep='\t')
## selected_metrics <- c("AT_DROPOUT","GC_DROPOUT","ON_BAIT_VS_SELECTED","PCT_PF_UQ_READS","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_50X")
## qcs.data <- qcs.data[which(col.names %in% selected_metrics)]

# CREATE A SCALED COPY OF THE DATA FRAME 
the.data.scaled <- pca.data
for (i in 2:ncol(the.data.scaled)) {
    mini <- min(the.data.scaled[,i])
    maxi <- max(the.data.scaled[,i])
    the.data.scaled[,i] <- apply(the.data.scaled, 1, function(row) { 
		row[[i]] <- (as.numeric(row[[i]]) - mini) / (maxi - mini)
	} )
}

# GET K-NEAREST NEIGHBORS FOR EACH SAMPLE
k.param <- 20
knns <- get.knn(the.data.scaled[,c(seq(2,ncol(the.data.scaled)))], k=k.param, algorithm="kd_tree")

# GENERATE A SINGLE FILE FOR EACH SAMPLE LISTING ITS K-NEAREST NEIGHBOR SAMPLE IDS
for (i in 1:nrow(the.data.scaled)) {
    fname <- paste(the.data.scaled$SAMPLE[i], ".", k.param, "nns.txt", sep="")
    nn.sampleids <- the.data.scaled$SAMPLE[ knns$nn.index[i,] ]
    write.table(nn.sampleids, fname, quote=F, row.names=F, col.names=F)
}

## system('for i in *.nns.txt; do sed \'s/$/.norm.cov.bed/\' $i > ${i%.norm.cov.bed}.ref.panel.files.txt; done')
system(paste("for i in *nns.txt; do sed \'s/$/.norm.cov.bed/\' $i | grep -f - ", in.sex, " > ${i%.20nns.txt}.ref.panel.files.txt; done"))

# TO CHECK HOW WELL EACH SAMPLE'S KNNS FIT, COMPUTE THE DISTANCE TO ITS KNN CLUSTER MEAN
the.data.scaled$DistanceToClusterMean <- sapply(1:nrow(the.data.scaled),  function(x) {
    this.knns <- knns$nn.index[x,];
    center <- colMeans(the.data.scaled[this.knns, 2:ncol(the.data.scaled)]);
    return(as.numeric(dist(rbind(as.numeric(the.data.scaled[x, 2:ncol(the.data.scaled)]), as.numeric(center)))))
})

## PLOT DISTANCE DISTRIBUTION
pdf("distance_distribution.pdf")
plot(ecdf(the.data.scaled$DistanceToClusterMean))
dev.off()
