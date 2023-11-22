library(ggtree)
library(dplyr)
library(Biostrings)
library(msa)
library(ape)
library(ips)
library(ggplot2)

fasta <- readDNAStringSet('~/Salba_RNA/genelists/mads_homologs/combined_dedup.fasta')
m <- msa(fasta, method = 'ClustalW')
mape <- as.DNAbin(m)
dist_mape <- dist.dna(mape,pairwise.deletion=TRUE)
dist_mape <- dist.dna(mape,model = "LOGDET")
dist_mape <- dist(mape)
any(is.na(dist_mape))  # Check for NAs
any(dist_mape ==0)    # Check for 0 values
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Paired"))(8)
t <- njs(dist_mape)
options(ignore.negative.edge=TRUE)
ggt <- ggtree(t, cex = 0.8, aes(color = branch.length)) + 
     scale_color_continuous(high = "#A6CEE3",low = "#f6b26b") +
     geom_tiplab(align = TRUE, size = 5) +
     geom_tippoint(size=3, fill="white", color="#e1ac0b", shape=21)+
     geom_treescale(y = -5, color = "#f1c232", fontsize = 10)
