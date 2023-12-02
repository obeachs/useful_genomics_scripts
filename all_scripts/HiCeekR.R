library(HiCeekR)
setwd("/Volumes/seagatedrive/")
HiCeekR()

install.packages(c("gplots", "randomizeBE", "MASS"))

install.packages("/Volumes/seagatedrive/HiCdatR_0.99.0.tar.gz", repos=NULL, type = "source")
library(HiCdatR)
f.source.organism.specific.code('/Volumes/seagatedrive/HiCdat/tair_10_HiCdat_R_code.R')
dataMatrix <- f.load.one.sample(dataDir = '/Volumes/seagatedrive/HiCdat/', files = c('read_mapped_to_fragments_matrix_with_counts_per_fragments.txt'),binSize = 1e6,repetitions = 50)
#dataMatrix_no_bin_size <- f.load.one.sample(dataDir = '/Volumes/seagatedrive/HiCdat/', files = c('read_mapped_to_fragments_matrix_with_counts_per_fragments.txt'),binSize = 0,repetitions = 50)

annotation <- f.read.annotation.via.fragment.annotation( annotationFile = "/Volumes/seagatedrive/HiCdat/newer_fragments_file.txt", binSize = 1e6,useLog = TRUE)
normalizedDataMatrix <- f.normalize.like.hu(dataMatrix = dataMatrix, binSize = 1e6, lenCol = 'length', gccCol = 'gcContent', mapCol = 'mappability', useNegativeBinomial = FALSE,annotation = annotationTable)
