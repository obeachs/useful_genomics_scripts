source('/Volumes/sesame/ALP_Omics/ChIP/validations/alp_visualisation_scripts.r')
col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/col-0_me3_annotated_midpeak_transposons_fixed_ref6_control_peaks_stats') %>% mutate(exp='Col-0')
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/alp2_me3_annotated_midpeak_transposons_fixed_ref6_control_peaks_stats') %>% mutate(exp='alp2-1')
ref6_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/ref6_alp2_me3_annotated_midpeak_transposons_fixed_ref6_control_peaks_stats') %>% mutate(exp='ref6_alp2-1')
ref6 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/stats/ref6_me3_annotated_midpeak_transposons_fixed_ref6_control_peaks_stats') %>% mutate(exp='ref6')

genes_of_interest <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps_with_clf_double/ref6_alp2_ref_control_overlap_with_rescues.csv')
#Usage example for plot_peaks
glist <- list(alp2, ref6, ref6_alp2, col)
titles <- list('alp2', 'ref6', 'ref6 alp2','Col-0')

for(i in 1:length(genes_of_interest$feature)){
  print(i)
  start <- genes_of_interest$start[i]
  end <- genes_of_interest$end[i]
  gene <- genes_of_interest$feature[i]
  print(typeof(start))
  print(typeof(end))
  print(typeof(gene))
  chr <- as.numeric(gsub('Chr','',genes_of_interest$chr[i]))
  print(typeof(chr))
 plot_peaks(glist,chr,start-2000,end + 2000,end, titles, paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps_with_clf_double/', gene,'.png'))
}