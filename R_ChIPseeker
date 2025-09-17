# Adapted from ChIPseeker https://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm39.refGene")
BiocManager::install("org.Mm.eg.db")

library(ChIPseeker) # version 1.40.0
library(TxDb.Mmusculus.UCSC.mm39.refGene) # version 3.19.0
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
library(org.Mm.eg.db) # version 3.19.1

#_________________________________________________________________________________
#  ** Schnetz 2010 HIGH thresh**
setwd("F:/")

peak2010 <- readPeakFile("GSM558674_CHD7_peak_High_thresh.txt")
covplot(peak2010)
peakAnno2010 <- annotatePeak(peak2010, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

plotAnnoPie(peakAnno2010)
plotDistToTSS(peakAnno2010)

symboldata2010 <- peakAnno2010@anno@elementMetadata@listData[["SYMBOL"]]
write.csv(symboldata2010, "Schnetz 2010__CHD7_peak_High_thresh.Anno.csv")
