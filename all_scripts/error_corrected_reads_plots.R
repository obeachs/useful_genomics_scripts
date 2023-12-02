library(seqinr)
library(adegenet)
library(ggplot2)
library(languageserver)
library(ggmsa)
library(ggpubr)
library(seqinr)
library(ggplot2)
library(patchwork)

nano <- read.delim('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_raw/all_nanopores_combind_1000bp.txt')
nano$Group <- 'Nanopore'
pac <- read.delim('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/5000bp_SRR1799170_used/joi_readlength_counts.txt')
pac$Group <- 'PacBio'

mix <- dplyr::full_join(nano,pac)




fa <- read.delim(fasta)
fa[fa == 0]  <- NA
fa <- na.omit(fa)

fa2 <- dplyr::slice(fa, -c(1,2,3,4,5,6,7,8,9))

p <- ggplot(data=mix, aes(x=Length, y=Count, fill = Group)) +
  geom_bar(stat="identity",position="dodge")
ggsave('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_raw/readlength_counts_over_1000bp_stack.png',p)  
p
