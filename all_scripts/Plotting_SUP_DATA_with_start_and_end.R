library(dplyr)
library(topGO)
library(KEGGREST)
library(org.At.tair.db)
library(GO.db)
library(ggplot2)
library(ggrepel)
SUP3.5D <- read.csv("/Volumes/seagatedrive/Work_folder/microarray_SUP/temp_SUP/Final_data/merge3.5D_sorted_0.05.txt")
SUP5D <- read.csv("/Volumes/seagatedrive/Work_folder/microarray_SUP/temp_SUP/Final_data/merge5D_sorted_0.05.txt")
load("~/Documents/annotation_TAIR10.Rdata")
setwd("/Volumes/seagatedrive/Work_folder/microarray_SUP/temp_SUP/Final_data")
  gff <- read.delim("/Volumes/seagatedrive/Work_folder/microarray_SUP/temp_SUP/TAIR10_GFF3_genes.txt")
#This is getting the entire list of genes that were differentially regulated in the SUP
#mutants
AT3_list <-SUP3.5D %>% full_join(SUP5D, by = "locus")
gff <- gff[,c(1,4,5,9)]
#Extracting only the TAIR IDs from the gff files - using a different package
#called stringr (needs to be called with stringr::)
TAIRIDs <- stringr::str_extract(gff$ID2, "AT.{1,7}")
gff$ID2 <- TAIRIDs
length(gff$ID2)
#Now we have to match the results from the 
matches <- gff$ID2 %in% AT3_list$locus
gff <- gff[matches,]
gff <- gff[match(unique(gff$ID2), gff$ID2),]
gff_AT3 <- gff[grep("AT3", gff$ID2),]

centromere_loc <- 13750000
names(AT3_list)[names(AT3_list) == "locus"] <- "ID2"
SUP_pos <- grep("AT3G23130", gff_AT3$ID2)
gff_AT3$midpoint <- (gff_AT3$END +gff_AT3$START)/2
gff_AT3$distancefromcentro <- gff_AT3$midpoint - centromere_loc
gff_AT3$allstartminusSUPstart <- gff_AT3$START - gff_AT3$START[SUP_pos]
gff_AT3$allendminusSUPstart <- gff_AT3$END - gff_AT3$START[SUP_pos]
gff_AT3$allstartminusSUPend <- gff_AT3$START - gff_AT3$END[SUP_pos]
gff_AT3$allendminusSUPend <- gff_AT3$END - gff_AT3$END[SUP_pos]
oopmerge <- merge(gff_AT3, AT3_list)
oopmerge <- oopmerge[,c(1,2,3,4,5,6,7,8,10)]
oopmerge <- oopmerge[,c(1,2,3,4,5,9)]
oopmerge <- oopmerge[-which(is.na(oopmerge$logFC.x)==TRUE),]
#Now going to compare the positions of the differentially expressed genes
#to SUPERMAN
#Adding genes at the start and the end of the chromosomes to show the 
#extent of the clustering
last_gene <- list(ID2="AT3G63460", Chromosome="Chr3", START=23430642, END=2347374,allstartminusSUPstart=15188340,logFC.x = 0)
first_gene <- list(ID2="AT3G01010", Chromosome="Chr3", START=4342, END=4818,allstartminusSUPstart=-8237960,logFC.x = 0)
oopmerge = rbind(oopmerge,last_gene, stringsAsFactors=FALSE)
oopmerge = rbind(oopmerge, first_gene, stringsAsFactors=FALSE)
ggplot(data = oopmerge, mapping = aes(x = distancefromcentro, y = logFC.x, color = ID2, label = ID2)) + 
  geom_point(alpha = 3) +
  geom_hline(yintercept =0) +
  geom_text_repel(size = 4, force = 2) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10)) +
  labs(x = "Distance from SUP locus", y = "log2 Fold Change")

oopmerge$ID2

please <- read.csv("~/Desktop/microarray_SUP/temp_SUP/merge3.5D_sorted_fix_29-01-20_0.05_test.txt")

