
'''Script that can be used to generate alignment plot from clustalo output.
Decided to make two different functions for DNA and protein for readability.
Sometimes might require some editing of the width and height of the graphs in the 
ggplot2::ggsave(outname_file,pop,width=N, height = N, limitsize = FALSE)
'''

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
      #Facet_msa defines how many bp per row there are
      #Char_width only effects the size of the font withinthe alignment blocks
      #it does not change the look or size of the seq_names
      pop <- ggmsa::ggmsa(file, border = NA, seq_name = T, char_width = 0.5,color = 'Zappo_AA',font = 'mono')+ggmsa::facet_msa(25) 
      ggplot2::ggsave(outname_file, pop, width = 20, height = 90, limitsize = FALSE) 
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
      pop <- ggmsa::ggmsa(file, border = NA, seq_name = T, char_width = 0.5,color = 'Shapely_NT',font = 'mono')+ggmsa::facet_msa(35) 
      ggplot2::ggsave(outname_file, pop, width = 20, height = 300, limitsize = FALSE) 
}
 
