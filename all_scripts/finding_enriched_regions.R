library(ggbio)
library(rtracklayer)
library(CSAR)
library(org.At.tair.db)
library(GenomeGraphs)
library(bPeaks)
library(ggplot2)
library(data.table)
library(UpSetR)

clf_28_alp1 <- '~/Documents/bigwigs_to_transfer/clf28_alp1-1_H3K27me3_summary_bamcompare.bw'
clf_28_alp2 <- '~/Documents/bigwigs_to_transfer/clf28_alp2-1_H3K27me3_summary_bamcompare.bw'
clf_28 <- '~/Documents/bigwigs_to_transfer/clf28_H3K27me3_summary_bamcompare.bw'
col_0 <- '~/Documents/bigwigs_to_transfer/col-0_H3K27me3_summary_bamcompare.bw'


clf_alp1 <- '~/Documents/alp_all_annotated_me3_peaks_epic2/2019_clf-28_alp1-1_annotated.csv'
clf_alp2 <- '~/Documents/alp_all_annotated_me3_peaks_epic2/2019_clf-28_alp2-1_annotated.csv'
alp1 <- '~/Documents/alp_all_annotated_me3_peaks_epic2/2020_alp1-1_annotated.csv'
alp2 <- '~/Documents/alp_all_annotated_me3_peaks_epic2/2020_alp2-1_annotated.csv'
alp2_chr <- '~/Documents/alp_all_annotated_me3_peaks_epic2/2021_alp2_annotated.csv'
col_0 <- '~/Documents/alp_all_annotated_me3_peaks_epic2/2019_Col-0_annotated.csv'
clf28 <- '~/Documents/alp_all_annotated_me3_peaks_epic2/2019_clf-28_annotated.csv'

bigwig_compare_difference <- function(control_bigwigcompare, experiment_bigwigcompare){
  cont <- import.bw(control_bigwigcompare)
  exp <- import.bw(experiment_bigwigcompare)
  cont <- Repitools::annoGR2DF(cont)
  exp <- Repitools::annoGR2DF(exp)
  merged_df <- left_join(cont, exp, by=c('chr','start'))
  merged_df <- na.omit(merged_df)
  merged_df <- merged_df[merged_df$score.x > 0,]
  merged_df <- merged_df[merged_df$score.y > 0,]
  merged_df$difference <- abs(merged_df$score.y) - abs(merged_df$score.x)
  
  mean(merged_df$difference)
  }

bigwig_compare_difference(clf_28,col_0)
bigwig_compare_difference(clf_28,clf_28_alp1)
bigwig_compare_difference(clf_28,clf_28_alp2)


#Checking for overlap with BR geneset published by Rahmani et al 2021
BR_genes <- '~/Documents/alp_all_annotated_me3_peaks_epic2/BR_genes.txt'

list.files('~/Documents/alp_all_annotated_me3_peaks_epic2/')

max_length <- length(genes)

force_df <- function(listy,maximum){
  c(listy,rep(NA,max_length - length(listy)))
}

peakfile <- read.csv(col_0)
genes <- as.data.frame(noquote(peakfile$feature))
geneslist <- read.csv(BR_genes, header = FALSE)

tpeakfile <- read.csv(clf_alp1)
tgenes <- as.data.frame(noquote(tpeakfile$feature))


col_0_int <- intersect(genes$`noquote(peakfile$feature)`, geneslist$V1)
alp1_int <- intersect(tgenes$`noquote(tpeakfile$feature)`,geneslist$V1)
intersect(col_0_int,alp1_int)

difg <- setdiff(genes,geneslist)

check_genes_present(BR_genes,clf28_alp1)
