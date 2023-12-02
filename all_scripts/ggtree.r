library(babette)
library(seqinr)
library(adegenet)
library(ggplot2)
library(ape)
library(languageserver)
library(ggmsa)
library(ggpubr)
library(seqinr)
library(ggplot2)
library(parallel)
library("KEGGREST")
library(Biostrings)
library(ggtree)
library(treeio)
library(phytools)

dna <- read.dna("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/phylogeny/AG_hits_same_lengths.mas", format = 'fasta', as.matrix = TRUE)
dist <- dist.gene(dna,)
nj <- nj(dist)
ggtree(nj)+ geom_tiplab()
ggsave('~/Desktop/tree_test_ATCG00120.1.pdf', ggtree(new) + geom_tiplab() + geom_tippoint(aes(color='red'), size=3, alpha=.75), width = 25)
dist
read.phylip("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/phylogeny/AG_hits_phylip.fa")

phylotools::read.phylip("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/phylogeny/AG_hits_phylip.fa")
new <- read.newick('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_corrected/phylogeny/ATCG00120.1.nwk')
ggtree(new)+ geom_tiplab()