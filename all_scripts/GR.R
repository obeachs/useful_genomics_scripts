library("GenomicRanges")
library("ChIPseeker")
library("ChIPpeakAnno")
library("chipenrich")
library ("rtracklayer")
library("TxDb.Athaliana.BioMart.plantsmart28")
library('ReactomePA')
txdb <- TxDb.Athaliana.BioMart.plantsmart28
#Making the peaks files
peakfile <- "/Volumes/seagatedrive/downloads/BT_endcomparison_narrowpeaks_edited_homered"
peak_df1 <- read.delim2(peakfile)
peaks_GR <-GRanges(seqnames = peak_df1[,'Chr'], IRanges(peak_df1[,'Start'],peak_df1[,'End']))
peakAnno <- annotatePeak(peaks_GR, tssRegion=c(-3000, 3000),TxDb=txdb)
plotAnnoPie(peakAnno)
enrichPathway(as.data.frame(peakAnno)$geneId)

#Setting up the ENSEMBL annoation files
txdb <- TxDb.Athaliana.BioMart.plantsmart28
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peaks_GR, windows = promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
