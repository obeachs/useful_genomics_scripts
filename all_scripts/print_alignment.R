
#Script that can be used to generate alignment plot from clustalo output.
#Decided to make two different functions for DNA and protein for readability.
#Sometimes might require some editing of the width and height of the graphs in the 
#ggplot2::ggsave(outname_file,pop,width=N, height = N, limitsize = FALSE)

library(dplyr)
library(RColorBrewer)
library(pals)
library(ggplot2)
protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
my_pal <- colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))
my_cutstom <- data.frame(names = c("A","C","-","G","T"), 
                         color = my_pal(5), 
                         stringsAsFactors = FALSE)
my_cutstom$color[3] <- '#ffffff'
color_df <- data.frame(
  stringsAsFactors = FALSE,
  names = c('*',"R", "H", "K", "D","E","S","T","N","Q","C","U","G","P","A","V","I","L","M","F","Y","W","-","B","J","Z","O"),
  color = c('#bdbdb5',"#d8f3cf", "#9ff381", "#009E73","#ebe94f","#faf88a","#FC9B7C","#FCA88B","#F34D37","#F75B40","#dcf3ff","#baf2ef","#a2d2df","#259aa1","#d2d4dc","#afafaf","#f8f8fa","#e5e6eb","#c0c2ce","#99a3ad","#879eb5","#c2c2c2","#ffffff","#ffffff","#ffffff","#ffffff","#ffffff")
)
color_df <- dplyr::left_join(my_cutstom, color_df,"names")


grn <- read.csv('~/Salba_RNA/genelists/trichome_GRN.csv')


AT4G18960_homologs_fixed_clustalo.fa
prot_align_to_pdf_new <- function(clustal_output) { 
      file <- clustal_output
      fasta <- seqinr::read.fasta(file)
      #All sequences in alignment fasta should be same length
      #Just taking the length of the first
      numseqs <- length(fasta)
      len <- length(fasta[[1]])
      outname_dir <- dirname(clustal_output)
      outname_file <- tools::file_path_sans_ext(clustal_output)
      outname_file <- paste(outname_file,'alignment_figure.pdf',sep = '_') 
      #outname_file <- '~/Desktop/testing_alignment.pdf'
      #Facet_msa defines how many bp per row there are
      #Char_width only effects the size of the font withinthe alignment blocks
      #it does not change the look or size of the seq_names
      pop <- ggmsa::ggmsa(file,seq_name = T, border =NA,char_width = 0.5,custom_color = color_df,font = 'mono')+ggmsa::facet_msa(25) 
      ggplot2::ggsave(outname_file, pop, width = 210, height = 297, unit='mm', limitsize = FALSE, dpi=500) 
}


dna_align_to_pdf_new <- function(clustal_output) { 
      file <- clustal_output
      fasta <- seqinr::read.fasta(file)
      #All sequences in alignment fasta should be same length
      #Just taking the length of the first
      numseqs <- length(fasta)
      len <- length(fasta[[1]])
      outname_dir <- dirname(clustal_output)
      outname_file <- tools::file_path_sans_ext(clustal_output)
      outname_file <- paste(outname_file,'alignment_figure.pdf',sep = '_') 
      #Facet_msa defines how many bp per row there are
      #Char_width only effects the size of the font withinthe alignment blocks
      #it does not change the look or size of the seq_names
      pop <- ggmsa::ggmsa(fasta, border = NA, seq_name = T,custom_color = my_cutstom, font = 'mono')+ggmsa::facet_msa(25) 
      ggplot2::ggsave(outname_file, pop, width = 210, height = 297, unit='mm', limitsize = FALSE, dpi=500) 
}

pop <-  ggmsa::ggmsa('/Users/josephbeegan/thesis_figs_and_tables/salba/homologs/prot/AT5G41315_pep_homologuesclustalo.fa',
seq_name = T, border ='grey',char_width = 0.5,custom_color = color_df,font = 'mono')+ggmsa::facet_msa(100) 

ggsave('/Users/josephbeegan/thesis_figs_and_tables/salba/homologs/prot/AT5G41315_pep_homologuesclustalo_alignment_figure.png', pop, width = 10,height=30,dpi=500, limitsize = F)
ggsave('~/thesis_figs_and_tables/salba/homologs/prot/AP1_homologs_pep.pdf', pop, width = 10, limitsize = F)





fasta <- seqinr::read.faswta('~/Salba_RNA/genelists/mads_homologs/homologs_fastas/AT2G46410_homologs_fixed_clustalo.fa')
pop <-  ggmsa::ggmsa('~/thesis_figs_and_tables/salba/homologs/cds/AT2G30432_cdna_homologs_clustalo.fa',
seq_name = T, border ='grey',char_width = 0.5,custom_color = my_cutstom,font = 'mono')+ggmsa::facet_msa(0) 

ggsave('~/thesis_figs_and_tables/salba/homologs/TCL1_homologs_pep.png', pop, width = 10,limitsize = F)
ggsave('~/thesis_figs_and_tables/salba/homologs/TCL1_homologs_pep.pdf', pop, width = 10, limitsize = F)     


files <- list.files('/Users/josephbeegan/Salba_RNA/gffcompare/novel_transcripts_fastas/', full.names = TRUE, pattern = 'dedup.fa')
names(files) <- file
for(file in files){
      print(file)
      for(i in grn$ID){
            if(grepl(i, file) != FALSE){
                  print(i)
            }
      }
}
dna_align_to_pdf_new('~/Salba_RNA/gffcompare/basic_gffcompare/novel_transcripts_fastas/AT2G30432_matches_dedup_clustalo.fa')
lapply(files, dna_align_to_pdf_new)
ggmsa::ggmsa(test,char_width = 0.5,custom_color = color_df,font = 'mono', seq_name = T)+ggmsa::facet_msa(25) 
prot_align_to_pdf_new('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/AT5G41315_pep_homologues_clustalo.fa')
prot_align_to_pdf_new("/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/final_figures/homologs/prot/AT5G41315_pep_homologuesclustalo.fa")
