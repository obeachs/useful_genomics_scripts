library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readAAStringSet("/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Gel_extracted_hom_mutants_and_wt_protein.fa", format = "fasta")

# look at some of the sequences (optional)
seqs
coolpatterns = c(alphabet(seqs, baseOnly=TRUE))
alpha
# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
#seqs <- OrientNucleotides(seqs)


# perform the alignment
aligned <- AlignSeqs(seqs)

writeXStringSet(aligned,
                file="/Volumes/seagatedrive/Project_folder/sinapis_assembly_shenanigans/alignments/TRY_sequences_aligned.fa")

# view the alignment in a browser (optional)
#For the colors, in order - GAPS,A,C,G,T

#IF DNA:
BrowseSeqs(seqs,colors=c("#CCCCCC",'#006080' ,"#FF6633", "#CC99FF", "#00CC66", "#0066CC",'#ff8080','#660000','#ffb84d','#995c00','#33ffd6','#00997a','#ff80d5','#99ff99','#006600','#ff4000','#ff6633','#e0b3ff','#7a00cc','#94b8b8','#4d94ff','#005ce6','#e6e600','#aaaa55','#002b80','#ffece6','#ff8c66','#d5ff80','#aaff00','#558000','#0099cc'),colWidth = 50,)
  
#IF PROTEINS
BrowseSeqs(seqs,colorPatterns = TRUE,patterns = coolpatterns, colors = substring(rainbow(length(coolpatterns),start=0.6, end=0.9,alpha = 0.2),1,7),colWidth = 60)
substrin
