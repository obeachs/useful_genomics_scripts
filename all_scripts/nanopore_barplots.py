library(seqinr)
library(adegenet)
library(ggplot2)
library(languageserver)
library(ggmsa)
library(ggpubr)
library(seqinr)
library(ggplot2)


file <- '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/joi_readlength_counts.txt'
fasta <- read.delim2(file)
fasta
p <- ggplot(data = fasta, aes(x=Length, y=Count)) + geom_bar(stat='identity')
ggsave('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/rjoi_readlength_counts.png',p)
max(fasta$Length)


