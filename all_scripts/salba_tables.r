library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gplots)
library(RVenn)
library(ggrepel)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(kableExtra)
library(tidyverse)
source('~/useful_genomics_scripts/all_scripts/formattable_functions_and_tables.R')


##QUAST results 
jgi <- read.csv('~/thesis_figs_and_tables/salba/JGI_quast.csv',header=T)
yang <- read.delim('~/thesis_figs_and_tables/salba/yang_quast.tsv', sep = '\t', header = T)
nanopore <- read.delim('~/thesis_figs_and_tables/salba/nanopore_SPAdes_assembly_quast.tsv', sep='\t', header = T)
tab <- rbind(jgi, yang, nanopore) %>% dplyr::select(X..contigs,Largest.contig,Total.length, 'GC....',N50, N90, L50,L90)
names(tab) <- c('Contigs','Largest Contig', 'Total length', 'GC Content', 'N50', 'N90','L50','L90') 
tab[] <- sapply(tab, as.numeric) 
tab <- as.data.frame(t(tab))
tab <- tab %>% mutate_all(~gsub("\\.00", "", as.character(.)))
colnames(tab)<- c('JGI Assembly', 'Yang et al 2023 Assembly', 'Nanopore reads assembly')

t <- kbl(tab, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down","add_linespace=-1.5mm")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:3, border_left = FALSE, border_right = FALSE)
cat(t)

##Alignment stats
stats <- read.csv('~/Salba_RNA/alignment_stats.csv') %>% mutate(Sample=gsub('_sorted','', Sample)) %>% mutate(Sample.1=gsub('_sorted','', Sample.1))
t <- kbl(stats, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1,3),bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:4, border_left = FALSE, border_right = FALSE)
cat(t)


##MADS homologs
species_counts <- read.csv('~/Salba_RNA/results/flower/mads/mads_homolog_counts.csv')
names(species_counts) <- gsub('_',' ', names(species_counts))
names(species_counts) <- gsub('\\.','-', names(species_counts))
t <- kbl(species_counts, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1,2),bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:4, border_left = FALSE, border_right = FALSE)
cat(t)


tair_mads <- read.delim('~/Salba_RNA/genelists/MADS.csv', header = F, sep='\t')
genes <- read.table('~/Salba_RNA/genelists/all_hits_strict.csv', header=T) %>% 
filter(tair %in% tair_mads$V2)
genes$gene_id <- gsub(",+$", "", genes$gene_id)
t <- kbl(genes, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1,2),bold = FALSE) %>% 
  kable_styling() %>%
  column_spec(1:4, border_left = FALSE, border_right = FALSE)
cat(t)
nnew_gennes <- genes %>% separate_rows('gene_id')
length(unique(nnew_gennes$gene_id))
