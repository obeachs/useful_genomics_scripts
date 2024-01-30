library(ggplot2)
library(dplyr)
library(ggVennDiagram)
library(gplots)
library(VennDiagram)
library(RVenn)
library(stringr)
library(ggpubr)
library(formattable)
library(ggrepel)
library(tidyverse)
library(ABuratinCircTarget)
library(circIMPACT)
library(kableExtra)
library(hypergea)
library(kableExtra)
library(gmp)
source('/Volumes/sesame/joerecovery/scripts/formattable_functions_and_tables.R')
source('/Volumes/sesame/ALP_Omics/ChIP/validations/alp_visualisation_scripts.r')
source('~/Salba_RNA/scripts/GO_dotplotter.R')
theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)

sample_table <- read.csv('~/Downloads/sample_ids.csv')
alp_data <- sample_table[!grepl('Input', sample_table$Condition),]
alp_data$Condition <- ifelse(grepl('H3', alp_data$Condition),alp_data$Condition, '') 
alp_data$ExperimentType <- paste(alp_data$Condition, alp_data$ExperimentType,sep = ' ')
alp_data <- alp_data %>% dplyr::select(Genotype, ExperimentType, Library) %>% 
dplyr::rename(Experiment=ExperimentType, 'Library type'= Library) %>%
dplyr::filter(!grepl('hdp', Genotype)) %>% distinct()
row.names(alp_data) <- NULL
t <- kbl(alp_data, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down","add_linespace=-1.5mm")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:3, border_left = FALSE, border_right = FALSE)
cat(t)



force_df <- function(listy,maximum){
  c(listy,rep('-',maximum - length(listy)))
}

get_unique_comparisons <- function(data) {
  unique_values <- unique(data$genotype)
  num_values <- length(unique_values)
  
  comparisons <- list()
  
  for (i in 1:(num_values - 1)) {
    for (j in (i + 1):num_values) {
      comparisons[[length(comparisons) + 1]] <- c(unique_values[i], unique_values[j])
    }
  }
  
  return(comparisons)
}


make_violin_plots_K27 <- function(df_of_fold_changes, outprefix, neg_num_y, no_neg_num_y){
  #Give a dataframe of genotype, FC_KO - no transformation of the negative data
  if (length(unique(df_of_fold_changes))>2){
    comparisons <- get_unique_comparisons(df_of_fold_changes)
    d <- ggviolin(df_of_fold_changes, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
        stat_summary(
          fun = "median",
          geom = "point",
          size = 3,
          color = "red"
        ) +
        stat_summary(
          fun = "mean",
          geom = "point",
          size = 3,
          color = "black"
        ) +
        theme +
        scale_y_continuous(
          breaks = seq(round(min(df_of_fold_changes$FC_KO))-1 , max(df_of_fold_changes$FC_KO), by = 1),
          limits = c(0, max(df_of_fold_changes$FC_KO) + no_neg_num_y)
        ) +
        theme(
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank(),
          legend.position = "none"
        ) +
        ylab('H3K27me3 Fold Change') +
        stat_compare_means(
          comparisons = comparisons,
          method = "t.test",  
          #label = "p.signif",
          #size = 4,
          #step.increase = 0.2
        ) + geom_hline(yintercept = 1, linetype = "dotted", color = "black")
    
    neg_df <- df_of_fold_changes %>% mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO))
    d_neg <- ggviolin(neg_df, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
        stat_summary(
          fun = "median",
          geom = "point",
          size = 3,
          color = "red"
        ) +
        stat_summary(
          fun = "mean",
          geom = "point",
          size = 3,
          color = "black"
        ) +
        theme +
        scale_y_continuous(
          breaks = seq(round(min(neg_df$FC_KO))-1 , max(neg_df$FC_KO),by = 1),
          limits = c(-(neg_num_y), max(neg_df$FC_KO) + neg_num_y)
        ) +
        theme(
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank(),
          legend.position = "none"
        ) +
        ylab('H3K27me3 Fold Change') +
        stat_compare_means(
          comparisons = comparisons,
          method = "t.test",  # You can use other methods such as "wilcox.test", "anova", etc.
          #label = "p.signif",
          #size = 4,
          #step.increase = 0.2
        )
    ggsave(paste0(outprefix,'_no_neg.pdf'),d, width=10, height=10)
    ggsave(paste0(outprefix,'.pdf'),d_neg, width=10, height=10)   
   }
  else{
    d <- ggviolin(df_of_fold_changes, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
        stat_summary(
          fun = "median",
          geom = "point",
          size = 3,
          color = "red"
        ) +
        stat_summary(
          fun = "mean",
          geom = "point",
          size = 3,
          color = "black"
        ) +
        theme +
        scale_y_continuous(
          breaks = seq(round(min(df_of_fold_changes$FC_KO))-1 , max(df_of_fold_changes$FC_KO), by = 1),
          limits = c(0, max(df_of_fold_changes$FC_KO) + 2.5)
        ) +
        theme(
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank(),
          legend.position = "none"
        ) +
        ylab('H3K27me3 Fold Change') +
        + geom_hline(yintercept = 1, linetype = "dotted", color = "black")
    neg_df <- df_of_fold_changes %>% mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO))
    d_neg <- ggviolin(neg_df, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
        stat_summary(
          fun = "median",
          geom = "point",
          size = 3,
          color = "red"
        ) +
        stat_summary(
          fun = "mean",
          geom = "point",
          size = 3,
          color = "black"
        ) +
        theme +
        scale_y_continuous(
          breaks = seq(round(min(neg_df$FC_KO))-1 , max(neg_df$FC_KO),by = 1),
          limits = c(-5, max(neg_df$FC_KO) + 5.5)
        ) +
        theme(
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank(),
          legend.position = "none"
        ) +
        ylab('H3K27me3 Fold Change')
    ggsave(paste0(outprefix,'_no_neg.pdf'),d, width=10, height=10)
    ggsave(paste0(outprefix,'.pdf'),d_neg, width=10, height=10)   
  }
}


get_TE_statistics <- function(peaks_file, # this is strictly required (no default)
                                            chr_col,
                                            start_col,
                                            end_col,
                                            header = T, # this has a default option, not strictly required
                                            sep = "\t"){

  peaks <- read.delim(peaks_file, header = header, sep = sep) %>%
    data.frame() %>%
    dplyr::rename("chr" = chr_col, "start" = start_col, "end" = end_col)
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(
    peaks,
    ignore.strand = T, # you sure we need this TRUE?
    seqnames.field = "chr", # can we get rid of this?
    start.field = "start",
    end.field = "end"
  )
  
  # Seems to be a lot of differences in the different TE databases for arabiopdsis
  # Some databases have TE_genes separate to TEs and don't include them.
  # This method is a brute force workaround that should solve all of the issues
  # Unfortunately doesn't allow inclusion of gene type in the initial output
  # Need to go from GRanges object, to dataframe and back to Granges in order to merge
  # the two datasets togethe


  #This is just for adding the type of gene/transposable element to the peaks
  # only multiplying it to avoid any awkwardness with the left_join later on
  all_genes <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/TAIR10_all_gene_ids.csv')
  all_genes <- all_genes[rep(seq_len(nrow(all_genes)), each = 20), ]
  #Importing the more substantial TE datasets

  trans <- read.delim('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_GFF_genes_transposons_slim.gff')
  trans <- trans[c(1,4,5,7,9)]
  names(trans) <- c('chr','start','end','strand','gene_id')
  trans$gene_id <- substr(trans$gene_id,1,9)
  rownames(trans) <- paste(trans$gene_id, rownames(trans),sep = '_')
  # Crucial that the keep.extra.columns=TRUE - otherwise the gene ids will not be
  # brought across making the annotation pointless
  gr_anno <- GenomicRanges::makeGRangesFromDataFrame(trans,keep.extra.columns=TRUE)
  
  annotated <- ChIPpeakAnno::annotatePeakInBatch(
    gr,
    AnnotationData = gr_anno,
    output = "nearestLocation",
    PeakLocForDistance = "middle")
  annotated <- Repitools::annoGR2DF(annotated)
  
  final_table <- dplyr::left_join(peaks, annotated, by='start',relationship='many-to-many')%>% distinct()
  final_table$feature <-  substr(final_table$feature,1,9)
  write.csv(final_table, '~/Desktop/test.csv', quote=F,row.names=F)
  print(nrow(final_table))
  all_TEs <- trans$gene_id[grepl('TE', trans$gene_id)]
  print(length(all_TEs))
  all_non_TEs <- trans$gene_id[!grepl('TE', trans$gene_id)]
  print(length(all_non_TEs))
  TEs <- nrow(final_table %>% filter(grepl('TE', feature)) %>% filter(FC_KO >1))
  non_TEs <- nrow(final_table %>% filter(!grepl('TE', feature)) %>% filter(FC_KO >1))

  print(TEs)
  print(non_TEs)
  # Total number of genes in the genome
  # total_genes_in_genome <- length(unique(final_table$gene_id))  # Replace with the actual total number of genes in your genome
  # expected_proportion_TE <- 0.21
  print(fisher.test(matrix(c(TEs,(nrow(all_TEs)-TEs),non_TEs,(nrow(all_non_TEs)-non_TEs)),nrow=2,ncol=2), alternative = 'greater'))
  
  print('1.2x')
  TEs <- nrow(final_table %>% filter(grepl('TE', feature)) %>% filter(FC_KO >1.2))
  non_TEs <- nrow(final_table %>% filter(!grepl('TE', feature)) %>% filter(FC_KO >1.2))

  print(TEs)
  print(non_TEs)
  # Total number of genes in the genome
  # total_genes_in_genome <- length(unique(final_table$gene_id))  # Replace with the actual total number of genes in your genome
  # expected_proportion_TE <- 0.21
  print(fisher.test(matrix(c(TEs,(nrow(all_TEs)-TEs),non_TEs,(nrow(all_non_TEs)-non_TEs)),nrow=2,ncol=2), alternative = 'greater'))
  print('1.5x')
  TEs <- nrow(final_table %>% filter(grepl('TE', feature)) %>% filter(FC_KO >1.5))
  non_TEs <- nrow(final_table %>% filter(!grepl('TE', feature)) %>% filter(FC_KO >1.5))

  print(TEs)
  print(non_TEs)
  # Total number of genes in the genome
  # total_genes_in_genome <- length(unique(final_table$gene_id))  # Replace with the actual total number of genes in your genome
  # expected_proportion_TE <- 0.21
  print(fisher.test(matrix(c(TEs,(nrow(all_TEs)-TEs),non_TEs,(0.21*nrow(all_non_TEs)-non_TEs)),nrow=2,ncol=2), alternative = 'greater'))
  print('2x')
  TEs <- nrow(final_table %>% filter(grepl('TE', feature)) %>% filter(FC_KO >2))
  non_TEs <- nrow(final_table %>% filter(!grepl('TE', feature)) %>% filter(FC_KO >2))

  print(TEs)
  print(non_TEs)
  print(fisher.test(matrix(c(TEs,(nrow(all_TEs)-TEs),non_TEs,(nrow(all_non_TEs)-non_TEs)),nrow=2,ncol=2), alternative = 'greater'))

}

tab <- matrix(c(0.21*25000,))
files <- list.files('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/', 
pattern='.txt', full.names=T)
for (i in files){
  if(grepl('diff', i)){
    print(i)
    get_TE_statistics(i,chr_col='Chromosome',start_col='Start',end_col = 'End')
  }
}

get_TE_statistics('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28.txt',chr_col='Chromosome',
start_col='Start',end_col = 'End')





raw <- list.files('/Volumes/sesame/ALP_Omics//ChIP/validations//stats//peak_number_stats/raw_peaks/', full.names = T)
raw_test <- read.csv(raw[1])
raw_df <- matrix(ncol = length(names(raw_test)), nrow = 0)
colnames(raw_df) <- names(raw_test)
for(i in raw){

    df <- read.csv(i)
    rownames(df) <- basename(i)
    raw_df <- rbind(raw_df, df)
}
write.csv(raw_df,'/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/raw_peaks/all_stats.csv', row.names=F, quote=F)



diff <- list.files('/Volumes/sesame/ALP_Omics//ChIP/validations//stats//peak_number_stats/diff_peaks',patter='m',full.names = T)
diff_test <- read.csv(diff[1])
diff_df <- matrix(ncol = length(names(diff_test)), nrow = 0)
colnames(diff_df) <- names(diff_test)
for(i in diff){
    df <- read.csv(i)
    rownames(df) <- basename(i)
    diff_df <- rbind(diff_df, df)
}


write.csv(diff_df,'/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/all_stats.csv', quote=F, row.names=F)

clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv')
clfalp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1)
clfalp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1)
overlaps_and_venn_output(genes=genelist, out_file_prefix = '~/thesis_figs_and_tables/alp/col_control_clf_doubles', plot_title = 'clf28 and clf28 double mutants - col-0 reference')
clfalp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.5| FC_KO < 0.75)
clfalp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.5 | FC_KO < 0.75)
overlaps_and_venn_output(genes=genelist, out_file_prefix = '~/thesis_figs_and_tables/alp/clf_control_clf_doubles', plot_title = 'clf28 double mutants - clf28 reference')

clf_alp1_clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.5| FC_KO < 0.75)
clf_alp2_clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.5| FC_KO < 0.75)
clf_alp1_rna <- read.csv("/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp1-1.csv")
clf_alp2_rna <- read.csv("/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp2-1.csv")

clf_alp1_list <- list(chipseq=clf_alp1_clf$feature,rnaseq= clf_alp1_rna$gene_id)
overlaps_and_venn_output(genes=clf_alp1_list, out_file_prefix = '~/thesis_figs_and_tables/alp/clf_control_clf_alp1_rna_chip', 
plot_title = 'clf28 alp1 ChIP-seq and RNA-seq - clf-28 reference',text_size=15,label_size=7.5)

clf_alp2_list <- list(chipseq=clf_alp2_clf$feature,rnaseq= clf_alp2_rna$gene_id)
overlaps_and_venn_output(genes=clf_alp2_list, out_file_prefix = '~/thesis_figs_and_tables/alp/clf_control_clf_alp2_rna_chip', 
plot_title = 'clf28 alp2 ChIP-seq and RNA-seq - clf-28 reference',text_size=15,label_size=7.5)

alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1)
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1)








swn_down_and_bound <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/swn_down_and_bound_genes.csv')

genelist <- list(clf_alp1=clfalp1$feature, clf_alp2=clfalp2$feature,swn=swn_down_and_bound$x)
overlaps_and_venn_output(genes=genelist, out_file_prefix = '~/thesis_figs_and_tables/alp/clf_alp_doubles_2x_UP_FC_swn_down_and_bound_overlaps_new', plot_title = 'clf alp double mutants and swn', text_size = 10)


diff_df_col <- diff_df %>% filter(grepl('2020', rownames(diff_df)) | grepl('col_control', rownames(diff_df))) %>%
rownames_to_column('Genotype') %>% mutate(Genotype=gsub("^.*peak_stats_(.*?)\\..*$", "\\1", Genotype)) %>% 
mutate(Genotype=gsub('_diff_',' ', Genotype)) %>% mutate(Genotype=gsub('_option','', Genotype)) %>% mutate(Genotype=gsub('_fixed',' ', Genotype))%>%
mutate(Genotype=gsub('_',' ', Genotype)) %>% mutate(Genotype=gsub('annotated',' ', Genotype)) %>% mutate(Genotype=gsub('annotations',' ', Genotype)) %>%
mutate(Genotype=gsub('both','(both)', Genotype))  %>% mutate(Genotype=gsub('midpeak','(nearest)', Genotype)) %>%
mutate("Histone Mark"=ifelse(grepl('me3', Genotype),'H3K27me3','H3K27me1')) %>% mutate(Genotype=gsub('me3','', Genotype)) %>%
mutate(Genotype=gsub('me1','', Genotype)) %>% filter(!grepl('transposons', Genotype))
names(diff_df_col) <- gsub('_',' ',names(diff_df_col))
names(diff_df_col) <- str_to_title(names(diff_df_col))

make_nice_table(df=diff_df_col, outname='~/thesis_figs_and_tables/alp/Differential_peaks_table.pdf', vwidth=1000, vheight=200)
f <- formattable(diff_df_col, 
            align =c(rep('c', ncol(diff_df_col))), list(Genotype = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))

f <- formattable(alp_data, 
            align =c(rep('c', ncol(alp_data))), list(Genotype = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))
export_formattable(f,'~/thesis_figs_and_tables/alp/experiment_list_BT.png')


f <- formattable(mtcars, 
            align =c(rep('c', ncol(alp_data))), list(cyl = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic")))) %>% as.htmlwidget() %>% 
  htmlwidgets::saveWidget(file="test.html")

f <- formattable(mtcars)
export_formattable(f,'~/Desktop/test_plot.pdf')



raw_df_col <- raw_df %>% rownames_to_column('Genotype') %>% mutate(Genotype=gsub("^.*peak_stats_(.*?)\\..*$", "\\1", Genotype)) %>% 
filter(!grepl('txt', Genotype)) %>%mutate(Genotype=gsub('_raw_',' ', Genotype)) %>% mutate(Genotype=gsub('_option','', Genotype)) %>% mutate(Genotype=gsub('_fixed',' ', Genotype))%>%
mutate(Genotype=gsub('_',' ', Genotype)) %>% mutate(Genotype=gsub('annotated',' ', Genotype)) %>% mutate(Genotype=gsub('annotations',' ', Genotype)) %>%
mutate(Genotype=gsub('both','(both)', Genotype))  %>% mutate(Genotype=gsub('midpeak','(nearest)', Genotype)) %>%
mutate("Histone Mark"=ifelse(grepl('me3', Genotype),'H3K27me3','H3K27me1')) %>% mutate(Genotype=gsub('me3','', Genotype)) %>%
mutate(Genotype=gsub('me1','', Genotype)) %>% mutate(Genotype=gsub('transposons','',Genotype)) %>% filter(grepl("\\(",Genotype))
names(raw_df_col) <- gsub('_',' ',names(raw_df_col))
names(raw_df_col) <- str_to_title(names(raw_df_col))
f <- formattable(raw_df_col, 
            align =c(rep('c', ncol(raw_df_col))), list(Genotype = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))
export_formattable(f,'~/thesis_figs_and_tables/alp/raw_peaks_table.png')


diff_df_col <- diff_df %>% filter(grepl('clf_control', rownames(diff_df))) %>%
rownames_to_column('Genotype') %>% mutate(Genotype=gsub("^.*peak_stats_(.*?)\\..*$", "\\1", Genotype)) %>% 
mutate(Genotype=gsub('_diff_',' ', Genotype)) %>% mutate(Genotype=gsub('_option','', Genotype)) %>% mutate(Genotype=gsub('_fixed',' ', Genotype))%>%
mutate(Genotype=gsub('_',' ', Genotype)) %>% mutate(Genotype=gsub('annotated',' ', Genotype)) %>% mutate(Genotype=gsub('annotations',' ', Genotype)) %>%
mutate(Genotype=gsub('both','(both)', Genotype))  %>% mutate(Genotype=gsub('midpeak','(nearest)', Genotype)) %>%
mutate("Histone Mark"=ifelse(grepl('me3', Genotype),'H3K27me3','H3K27me1')) %>% mutate(Genotype=gsub('me3','', Genotype)) %>%
mutate(Genotype=gsub('me1','', Genotype)) %>% filter(!grepl('transposons', Genotype)) %>% mutate(Genotype=gsub('v clf28','', Genotype))
names(diff_df_col) <- gsub('_',' ',names(diff_df_col))
names(diff_df_col) <- str_to_title(names(diff_df_col))
f <- formattable(diff_df_col, 
            align =c(rep('c', ncol(diff_df_col))), list(Genotype = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))
export_formattable(f,'~/thesis_figs_and_tables/alp/Differential_peaks_clf_control_table.png')


diff_df_col <- diff_df %>% filter(grepl('ref_control', rownames(diff_df))) %>%
rownames_to_column('Genotype') %>% mutate(Genotype=gsub("^.*peak_stats_(.*?)\\..*$", "\\1", Genotype)) %>% 
mutate(Genotype=gsub('_diff_',' ', Genotype)) %>% mutate(Genotype=gsub('_option','', Genotype)) %>% mutate(Genotype=gsub('_fixed',' ', Genotype))%>%
mutate(Genotype=gsub('_',' ', Genotype)) %>% mutate(Genotype=gsub('annotated',' ', Genotype)) %>% mutate(Genotype=gsub('annotations',' ', Genotype)) %>%
mutate(Genotype=gsub('both','(both)', Genotype))  %>% mutate(Genotype=gsub('midpeak','(nearest)', Genotype)) %>%
mutate("Histone Mark"=ifelse(grepl('me3', Genotype),'H3K27me3','H3K27me1')) %>% mutate(Genotype=gsub('me3','', Genotype)) %>%
mutate(Genotype=gsub('me1','', Genotype)) %>% filter(!grepl('transposons', Genotype))%>% mutate(Genotype=gsub('ref6 control','', Genotype))
names(diff_df_col) <- gsub('_',' ',names(diff_df_col))
names(diff_df_col) <- str_to_title(names(diff_df_col))
f <- formattable(diff_df_col, 
            align =c(rep('c', ncol(diff_df_col))), list(Genotype = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))
export_formattable(f,'~/thesis_figs_and_tables/alp/Differential_peaks_ref_control_table.png')




shu <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv') %>% filter(Type.of.K27.reduction.pattern=='III')
bt <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO >2 | FC_KO < 0.5)
clf_up <-length(shu$clf_K27.increase[grepl('AT',shu$clf_K27.increase)])
clf_down <- length(shu$clf_K27.decrease[grepl('AT',shu$clf_K27.decrease)])
swn_up <- length(shu$swn_K27.increase[grepl('AT',shu$swn_K27.increase)])
swn_down <- length(shu$swn_K27.decrease[grepl('AT',shu$swn_K27.decrease)])
bt_up <- length(bt$feature[bt$FC_KO>1])
bt_down <- length(bt$feature[bt$FC_KO<1])


#Table of genes with changes in SWN H3k27me3
gene_names <- read.csv('~/Salba_RNA/genelists//all_hits_names.csv') %>% dplyr::select(tair, gene.name) %>% distinct()
swn_up <- as.data.frame(shu$swn_K27.increase[grepl('AT',shu$swn_K27.increase)]) 
names(swn_up) <- c('tair')
swn_up <- swn_up %>% left_join(gene_names) %>% mutate(gene.name=ifelse(is.na(gene.name), tair, gene.name)) %>% dplyr::rename('TAIR ID'=tair, 'Gene Name'=gene.name)
make_nice_table(swn_up, outname='~/thesis_figs_and_tables/alp/swn_up_gene_ids_and_names.png')

swn_down <- as.data.frame(shu$swn_K27.decrease[grepl('AT',shu$swn_K27.decrease)]) 
names(swn_down) <- c('tair')
swn_down <- swn_down %>% left_join(gene_names) %>% mutate(gene.name=ifelse(is.na(gene.name), tair, gene.name)) %>% dplyr::rename('TAIR ID'=tair, 'Gene Name'=gene.name)
make_nice_table(swn_down, outname='~/thesis_figs_and_tables/alp/swn_down_gene_ids_and_names.png')

clf_up <-as.data.frame(shu$clf_K27.increase[grepl('AT',shu$clf_K27.increase)])
names(clf_up) <- c('tair')
clf_up <- clf_up %>% left_join(gene_names) %>% mutate(gene.name=ifelse(is.na(gene.name), tair, gene.name)) %>% dplyr::rename('TAIR ID'=tair, 'Gene Name'=gene.name)
make_nice_table(clf_up, outname='~/thesis_figs_and_tables/alp/shu_clf_up_gene_ids_and_names.png')

clf_down <-as.data.frame(shu$clf_K27.decrease[grepl('AT',shu$clf_K27.decrease)])
names(clf_down) <- c('tair')
clf_down <- clf_down %>% left_join(gene_names) %>% mutate(gene.name=ifelse(is.na(gene.name), tair, gene.name)) %>% dplyr::rename('TAIR ID'=tair, 'Gene Name'=gene.name)
make_nice_table(clf_down[200:215,], outname='~/thesis_figs_and_tables/alp/shu_clf_down_gene_ids_and_names_215.png')


f <- formattable(df, 
            align =c(rep('c', ncol(df))), list(Experiment = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))
export_formattable(f,'~/thesis_figs_and_tables/alp/shu_bt_overall_H3K27me3.png')

df <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/H3K27me3_counts.csv')
names(df) <- c("Experiment","H3K27me3 UP", "H3K27me3 Down")

















###Prepping for presentation
clf_volcano_all <- clf %>% dplyr::select(FC_KO, FDR_KO) %>% 
 mutate(genotype=paste('clf28 all peaks (n =', nrow(clf_volcano_all),')'))
clf_volcano_strict <- clf %>% dplyr::select(FC_KO, FDR_KO) %>% filter(FC_KO >2 | FC_KO < 0.5)%>% 
mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>% mutate(genotype=paste('clf28 2x FC peaks (n =', nrow(clf_volcano_strict),')'))
clf_volcano <- rbind(clf_volcano_all, clf_volcano_strict)
tab <-  read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/2019_col_control_peak_stats_clf28_me3_diff_annotated_midpeak.csv')

plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/2019_col_control_peak_stats_clf28_me3_diff_annotated_midpeak.csv',
'~/thesis_figs_and_tables/alp/fold_change_barplots/clf28_col_control_H3K27me3_fold_changes_barplot.png')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/2020_peak_stats_alp1_me3_diff_annotated_midpeak.csv',
'~/thesis_figs_and_tables/alp/fold_change_barplots/alp1_col_control_H3K27me3_fold_changes_barplot.png')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/2020_peak_stats_alp2_me3_diff_annotated_midpeak.csv',
'~/thesis_figs_and_tables/alp/fold_change_barplots/alp2_col_control_H3K27me3_fold_changes_barplot.png')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats//diff_peaks/2019_col_control_peak_stats_clf28_alp1_me3_diff_annotated_both_option.csv',
'~/thesis_figs_and_tables/alp/fold_change_barplots/clf_alp1_col_control_H3K27me3_fold_changes_barplot.png')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats//diff_peaks/2019_col_control_peak_stats_clf28_alp2_me3_diff_annotated_both_option.csv',
'~/thesis_figs_and_tables/alp/fold_change_barplots/clf_alp2_col_control_H3K27me3_fold_changes_barplot.png')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_clf_v_col_merged_bams/peak_stats_clf_merged_me3-and-col_merged_me3-W200-G600-summary_unedited',
'~/thesis_figs_and_tables/alp/fold_change_barplots/clf29_col_control_H3K27me3_fold_changes_barplot.png')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_swn_v_col_merged_bams/peak_stats_swn_merged_me3-and-col_merged_me3-W200-G600-summary_unedited',
'~/thesis_figs_and_tables/alp/fold_change_barplots/swn4_col_control_H3K27me3_fold_changes_barplot.png')

ggsave('~/thesis_figs_and_tables/alp/clf28_col_control_H3K27me3_fold_changes_volcanoplot_no_neg.pdf',d, width=10, height=10)
ggsave('~/thesis_figs_and_tables/alp//fold_change_barplots/clf28_col_control_H3K27me3_fold_changes_barplot.pdf',p, width=10, height=10)


alp_1_single <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_midpeak.csv') %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(genotype=paste0('alp1 all peaks (n = ', nrow(alp_1_single),')'))
alp_1_strict<- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO >2 | FC_KO < 0.5) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>% mutate(genotype=paste0('alp1 2x FC (n = ', nrow(alp_1_strict),')'))
alp1_volcano <- rbind(alp_1_single, alp_2_single) 
alp1_volcano <- alp1_volcano %>%  filter(FC_KO < 15 & FC_KO > -15)
d <- ggviolin(alp1_volcano, x = "genotype", y = "FC_KO",palette = "Paired", # Add jittered points for better visualization
  fill = "genotype")+    # Fill the violins with colors based on the "exp" variable)
  stat_summary(
    fun = "median", # Calculate median for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of median points
    color = "red"   # Set the color of median points
  ) +
   theme+
   scale_y_continuous(breaks = seq(round(min(alp1_volcano$FC_KO))-1 , max(alp1_volcano$FC_KO) -10,by = 1),
    limits = c(-0, max(alp1_volcano$FC_KO) -10))+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),legend.position = "none"
  ) + ylab('H3K27me3 Fold Change Relative to Col-0')+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")
tab <-  read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/2020_peak_stats_alp1_me3_diff_annotated_midpeak.csv')
data_long <- gather(tab, key, value)%>% filter(!grepl('total', key))

# Extract the direction (up or down) from the key column
data_long$direction <- ifelse(grepl("up", data_long$key), "H3K27me3 Up", "H3K27me3 Down")
data_long$key <- gsub('peaks_','', data_long$key)
data_long$key <- gsub('up_','', data_long$key)
data_long$key <- gsub('down_','', data_long$key)
data_long$key <- gsub('up','All', data_long$key)
data_long$key <- gsub('down','All', data_long$key)
data_long$key <- gsub(' ','', data_long$key)
data_long$key <- stringr::str_to_title(data_long$key)

p <- ggplot(data_long, aes(x = key, y = value, fill = direction, label = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("H3K27me3 Down" = "#E1BE6A", "H3K27me3 Up" = "#40B0A6")) +
  theme_minimal() +
  xlab('Peak Fold Change relative to Col-0') +
  ylab('Count') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 20),
    axis.ticks.x = element_blank(),  # Remove ticks
    panel.grid.major.x = element_blank(),  # Remove gridlines
    panel.grid.minor.x = element_blank(),  # Remove gridlines
    legend.title = element_text(size = 15),  # Change legend title and size
    legend.text = element_text(size = 15)
  ) +
  geom_text(
    size = 5,
    position = position_dodge(width = 0.9),  # Adjust width as needed
    vjust = -0.5,  # Adjust vjust as needed
    colour = "black"
  ) +theme+
  scale_x_discrete(limits=c('All', '1.2x', '1.5x','2x','3x'))+
  labs(fill = "Direction of H3K27me3 Change")
  
ggsave('~/thesis_figs_and_tables/alp/alp1_col_control_H3K27me3_fold_changes_volcanoplot_no_neg.pdf',d, width=10, height=10)
ggsave('~/thesis_figs_and_tables/alp/alp1_col_control_H3K27me3_fold_changes_barplot.pdf',p, width=10, height=10)

alp_2_single <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_midpeak.csv') %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(genotype=paste0('alp2 all peaks (n = ', nrow(alp_2_single),')'))
alp_2_strict <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO >2 | FC_KO < 0.5) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>% mutate(genotype=paste0('alp2 2x FC (n =', nrow(alp_2_strict),')'))
alp2_volcano <- rbind(alp_2_single, alp_2_strict)
alp2_volcano <- alp2_volcano %>%  filter(FC_KO < 15 & FC_KO > -15)
d <- ggviolin(alp2_volcano, x = "genotype", y = "FC_KO",palette = "Paired", # Add jittered points for better visualization
  fill = "genotype")+    # Fill the violins with colors based on the "exp" variable)
  stat_summary(
    fun = "median", # Calculate median for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of median points
    color = "red"   # Set the color of median points
     ) +
   theme+
   scale_y_continuous(breaks = seq(round(min(alp2_volcano$FC_KO))-1, round(max(alp2_volcano$FC_KO))+1, by = 1))+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),legend.position = "none"

  ) + ylab('H3K27me3 Fold Change Relative to Col-0')
tab <-  read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/2020_peak_stats_alp2_me3_diff_annotated_midpeak.csv')
data_long <- gather(tab, key, value)%>% filter(!grepl('total', key))

# Extract the direction (up or down) from the key column
data_long$direction <- ifelse(grepl("up", data_long$key), "H3K27me3 Up", "H3K27me3 Down")
data_long$key <- gsub('peaks_','', data_long$key)
data_long$key <- gsub('up_','', data_long$key)
data_long$key <- gsub('down_','', data_long$key)
data_long$key <- gsub('up','All', data_long$key)
data_long$key <- gsub('down','All', data_long$key)
data_long$key <- gsub(' ','', data_long$key)
data_long$key <- stringr::str_to_title(data_long$key)

p <- ggplot(data_long, aes(x = key, y = value, fill = direction, label = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("H3K27me3 Down" = "#E1BE6A", "H3K27me3 Up" = "#40B0A6")) +
  theme_minimal() +
  xlab('Peak Fold Change relative to Col-0') +
  ylab('Count') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 20),
    axis.ticks.x = element_blank(),  # Remove ticks
    panel.grid.major.x = element_blank(),  # Remove gridlines
    panel.grid.minor.x = element_blank(),  # Remove gridlines
    legend.title = element_text(size = 15),  # Change legend title and size
    legend.text = element_text(size = 15)
  ) +
  geom_text(
    size = 5,
    position = position_dodge(width = 0.9),  # Adjust width as needed
    vjust = -0.5,  # Adjust vjust as needed
    colour = "black"
  ) +x+
  scale_x_discrete(limits=c('All', '1.2x', '1.5x','2x','3x'))+
  labs(fill = "Direction of H3K27me3 Change")

ggsave('~/thesis_figs_and_tables/alp/alp2_col_control_H3K27me3_fold_changes_volcanoplot.pdf',d, width=10, height=10)
ggsave('~/thesis_figs_and_tables/alp/alp2_col_control_H3K27me3_fold_changes_barplot.pdf',p, width=10, height=10)
write.csv(clf_alp2,'~/Desktop/clf_alp2.csv', quote=F, row.names=F)

names <- read.csv('~/Salba_RNA/genelists/all_hits_names.csv') %>% dplyr::select(tair, gene.name) %>% distinct() %>% dplyr::rename(feature=tair)
clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv')
clf_down <- clf %>% filter(FC_KO < 1)
clf_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_midpeak.csv')
clf_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv')


clf_alp1_clf_down <- clf_alp1 %>% filter(feature %in% clf_down$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>% mutate(genotype=paste0('clf28 alpl1 FC (n = ', nrow(clf_alp1_clf_down),')')) %>%
mutate(exp='clf_alpl1')
clf_alp2_clf_down <- clf_alp2 %>% filter(feature %in% clf_down$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% 
 mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>% mutate(genotype=paste0('clf28 alpl2 FC (n = ', nrow(clf_alp2_clf_down),')')) %>%
mutate(exp='clf_alpl2')
clf_down <- clf_down %>% filter(FC_KO < 1) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>% mutate(genotype=paste0('clf28 FC (n = ', nrow(clf_down),')')) %>%
mutate(exp='clf')

clf_alp1_clf_down <- clf_alp1 %>% filter(feature %in% clf_down$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(genotype=paste0('clf28 alpl1 FC (n = ', nrow(clf_alp1_clf_down),')')) %>%
mutate(exp='clf_alpl1')
clf_alp2_clf_down <- clf_alp2 %>% filter(feature %in% clf_down$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(genotype=paste0('clf28 alpl2 FC (n = ', nrow(clf_alp2_clf_down),')')) %>%
mutate(exp='clf_alpl2')
clf_down <- clf_down %>% filter(FC_KO < 1) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(genotype=paste0('clf28 FC (n = ', nrow(clf_down),')')) %>%
mutate(exp='clf')

doubles_down <- rbind(clf_down,clf_alp1_clf_down,clf_alp2_clf_down )

my_comparisons <- list(c(unique(doubles_down$genotype)[1],unique(doubles_down$genotype)[2]),
c(unique(doubles_down$genotype)[2],unique(doubles_down$genotype)[3]),c(unique(doubles_down$genotype)[1],unique(doubles_down$genotype)[3]))

d <- ggviolin(doubles_down, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme +
  scale_y_continuous(
    breaks = seq(round(min(doubles_down$FC_KO)) - 1, round(max(doubles_down$FC_KO)) + 1, by = 1),
    limits = c(0, max(doubles_down$FC_KO) + 2)
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab('H3K27me3 Fold Change Relative to Col-0') +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",  # You can use other methods such as "wilcox.test", "anova", etc.
    #label = "p.signif",
    #size = 4,
    #step.increase = 0.2
  ) + geom_hline(yintercept = 1, linetype = "dotted", color = "black")

ggsave('~/thesis_figs_and_tables/alp/clf_alp_double_clf_H3K27me3_down_shared_peaks_violin_no_neg.pdf',d, height=10, width = 10)


clf_up <- clf %>% filter(FC_KO > 1)
clf_alp1_clf_up <- clf_alp1 %>% filter(feature %in% clf_up$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>% mutate(genotype=paste0('clf28 alpl1 2x FC (n = ', nrow(clf_alp1_clf_up),')')) %>%
mutate(exp='clf_alpl1')
clf_alp2_clf_up <- clf_alp2 %>% filter(feature %in% clf_up$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>%  mutate(genotype=paste0('clf28 alpl2 2x FC (n = ', nrow(clf_alp2_clf_up),')')) %>%
mutate(exp='clf_alpl2')
clf_up <- clf_up %>% filter(FC_KO >1) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(FC_KO = ifelse((FC_KO < 1),-1*(1/FC_KO),FC_KO)) %>% mutate(genotype=paste0('clf28 2x FC (n = ', nrow(clf_up),')')) %>%
mutate(exp='clf')
doubles_up <- rbind(clf_up,clf_alp1_clf_up,clf_alp2_clf_up )

my_comparisons <- list(c(unique(doubles_up$genotype)[1],unique(doubles_up$genotype)[2]),
c(unique(doubles_up$genotype)[2],unique(doubles_up$genotype)[3]),c(unique(doubles_up$genotype)[1],unique(doubles_up$genotype)[3]))

make_violin_plots_K27(doubles_up,'~/Desktop/test_doubles')


d <- ggviolin(doubles_up, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme +
  scale_y_continuous(
    breaks = seq(round(min(doubles_up$FC_KO))-1 , max(doubles_up$FC_KO), by = 1),
    limits = c(0, max(doubles_up$FC_KO) + 2.5)
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab('H3K27me3 Fold Change Relative to Col-0') +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",  # You can use other methods such as "wilcox.test", "anova", etc.
    #label = "p.signif",
    #size = 4,
    #step.increase = 0.2
  ) + geom_hline(yintercept = 1, linetype = "dotted", color = "black")

d_neg <- ggviolin(doubles_up, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme +
  scale_y_continuous(
    breaks = seq(round(min(doubles_up$FC_KO))-1 , max(doubles_up$FC_KO),by = 1),
    limits = c(-5, max(doubles_up$FC_KO) + 5.5)
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab('H3K27me3 Fold Change Relative to Col-0') +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",  # You can use other methods such as "wilcox.test", "anova", etc.
    #label = "p.signif",
    #size = 4,
    #step.increase = 0.2
  )

ggsave('~/thesis_figs_and_tables/alp/clf_alp_double_clf_H3K27me3_UP_FC_shared_peaks_violin.pdf',d_neg, height=10, width = 10)



clf_all <- clf
clf_alp1_clf_all <- clf_alp1 %>% filter(feature %in% clf_all$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(genotype=paste0('clf28 alpl1 (n = ', nrow(clf_alp1_clf_up),')')) %>%
mutate(exp='clf_alpl1')
clf_alp2_clf_all <- clf_alp2 %>% filter(feature %in% clf_all$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(genotype=paste0('clf28 alpl2(n = ', nrow(clf_alp2_clf_up),')')) %>%
mutate(exp='clf_alpl2')
clf_all <- clf_all %>% dplyr::select(FC_KO, FDR_KO) %>% 
mutate(genotype=paste0('clf28(n = ', nrow(clf_up),')')) %>%
mutate(exp='clf')
doubles_all <- rbind(clf_all,clf_alp1_clf_all,clf_alp2_clf_all)

my_comparisons <- list(c(unique(doubles_all$genotype)[1],unique(doubles_all$genotype)[2]),
c(unique(doubles_all$genotype)[2],unique(doubles_all$genotype)[3]),c(unique(doubles_all$genotype)[1],unique(doubles_all$genotype)[3]))

d <- ggviolin(doubles_all, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme +
  scale_y_continuous(
    breaks = seq(round(min(doubles_all$FC_KO))-1 , max(doubles_all$FC_KO), by = 1),
    limits = c(0, max(doubles_all$FC_KO) + 2.5)
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab('H3K27me3 Fold Change Relative to Col-0') +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",  # You can use other methods such as "wilcox.test", "anova", etc.
    #label = "p.signif",
    #size = 4,
    #step.increase = 0.2
  ) + geom_hline(yintercept = 1, linetype = "dotted", color = "black")

d_neg <- ggviolin(doubles_all, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme +
  scale_y_continuous(
    breaks = seq(round(min(doubles_all$FC_KO))-1 , max(doubles_all$FC_KO),by = 1),
    limits = c(-5, max(doubles_all$FC_KO) + 5.5)
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab('H3K27me3 Fold Change Relative to Col-0') +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",  # You can use other methods such as "wilcox.test", "anova", etc.
    #label = "p.signif",
    #size = 4,
    #step.increase = 0.2
  )

ggsave('~/thesis_figs_and_tables/alp/clf_alp_double_clf_H3K27me3_all_shared_peaks_violin_no_negs.pdf',d, height=10, width = 10)




clf_alp1_clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv') %>%
dplyr::select(FC_KO, FDR_KO)  %>% 
mutate(genotype=paste0('clf28 alpl1 (n = ', nrow(clf_alp1_clf),')')) %>%mutate(exp='clf_alpl1')
clf_alp2_clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv') %>%
dplyr::select(FC_KO, FDR_KO)  %>% 
mutate(genotype=paste0('clf28 alpl2  (n = ', nrow(clf_alp2_clf),')')) %>% mutate(exp='clf_alpl2')
clf_control <- rbind(clf_alp1_clf,clf_alp2_clf)
my_comparisons <- list(c(unique(clf_control$genotype)[1],unique(clf_control$genotype)[2]))
d <- ggviolin(clf_control, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme + 
  scale_y_continuous(
    breaks = seq(round(min(clf_control$FC_KO))-1 , max(clf_control$FC_KO), by = 1),
    limits = c(0, max(clf_control$FC_KO) + 2.5)
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab('H3K27me3 Fold Change Relative to clf28') +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",  # You can use other methods such as "wilcox.test", "anova", etc.
    #label = "p.signif",
    #size = 4,
    #step.increase = 0.2
  ) + geom_hline(yintercept = 1, linetype = "dotted", color = "black")
d_neg <- ggviolin(clf_control, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
  stat_summary(
    fun = "median",
    geom = "point",
    size = 3,
    color = "red"
  ) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 3,
    color = "black"
  ) +
  theme +
  scale_y_continuous(
    breaks = seq(round(min(clf_control$FC_KO))-1 , max(clf_control$FC_KO),by = 1),
    limits = c(-5, max(clf_control$FC_KO)+1)
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab('H3K27me3 Fold Change Relative to clf28') +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",  # You can use other methods such as "wilcox.test", "anova", etc.
    #label = "p.signif",
    #size = 4,
    #step.increase = 0.2
  )

ggsave('~/thesis_figs_and_tables/alp/clf_alp_double_clf_control_violin.pdf',d_neg, height=10, width = 10)
ggsave('~/thesis_figs_and_tables/alp/clf_alp_double_clf_control_violin_no_neg.pdf',d, height=10, width = 10)

tab <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats//diff_peaks/2019_clf_control_peak_stats_clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats//diff_peaks/2019_col_control_peak_stats_clf28_alp2_me3_diff_annotated_both_option.csv',
'~/thesis_figs_and_tables/alp/fold_change_barplots/clf_alp2_col_control_H3K27me3_fold_changes_barplot.pdf')
ref6_alp2_stats <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/FC_stats', full.names=T)
plot_peak_stats_barplot <- function(stats_table, outname, control='Col-0'){

  tab <- read.csv(stats_table)
  data_long <- gather(tab, key, value)%>% filter(!grepl('total', key))
    
  # Extract the direction (up or down) from the key column
  data_long$direction <- ifelse(grepl("up", data_long$key), "H3K27me3 Up", "H3K27me3 Down")
  data_long$key <- gsub('peaks_','', data_long$key)
  data_long$key <- gsub('up_','', data_long$key)
  data_long$key <- gsub('down_','', data_long$key)
  data_long$key <- gsub('up','All', data_long$key)
  data_long$key <- gsub('down','All', data_long$key)
  data_long$key <- gsub(' ','', data_long$key)
  data_long$key <- stringr::str_to_title(data_long$key)

  p <- ggplot(data_long, aes(x = key, y = value, fill = direction, label = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("H3K27me3 Down" = "#E1BE6A", "H3K27me3 Up" = "#40B0A6")) +
    theme_minimal() +
    xlab(paste0('Peak Fold Change relative to ', control)) +
    ylab('Count') +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 20),
      axis.ticks.x = element_blank(),  # Remove ticks
      panel.grid.major.x = element_blank(),  # Remove gridlines
      panel.grid.minor.x = element_blank(),  # Remove gridlines
      legend.title = element_text(size = 15),  # Change legend title and size
      legend.text = element_text(size = 15)
    ) +
    geom_text(
      size = 5,
      position = position_dodge(width = 0.9),  # Adjust width as needed
      vjust = -0.5,  # Adjust vjust as needed
      colour = "black"
    ) +theme+
    scale_x_discrete(limits=c('All', '1.2x', '1.5x','2x','3x'))+
    labs(fill = "Direction of H3K27me3 Change")


  ggsave(outname,p, width=10, height=10, dpi=500)
}

for(i in 1:length(ref6_alp2_stats)){
  out <- gsub('.csv','.pdf', ref6_alp2_stats[i])
  if(grepl('ref6_control', ref6_alp2_stats[i])){
  plot_peak_stats_barplot(ref6_alp2_stats[i], out, control='ref6')
  }
  else{
    plot_peak_stats_barplot(ref6_alp2_stats[i], out)
  }
}

plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/sicer_ref6_alp2_v_UDIcol_merged_bams', 
'~/thesis_figs_and_tables/alp/fold_change_barplots/alp2_FC_ref_control_table.csvpeak_stats_alp2_me3_diff_ref6_control_annotated_midpeak.pdf', control='ref6')

##Donut plot of clf alp double mutant v clf H3K27me3
#getting clf_h3K27me3 down
clf_down <- clf %>% filter(FC_KO < 1)
clf_alp1_clf_down <- clf_alp1_clf %>% filter(feature %in% clf_down$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% filter(FC_KO >1)
df <- data.frame(
  Status = c("H3K27me3 Rescued", "H3K27me3 Not Rescued"),
  Count = c(nrow(clf_alp1_clf_down), nrow(clf_down)-nrow(clf_alp1_clf_down))
)
df$ratio <- round((df$Count/sum(df$Count))*100)
p <- donut(df)
ggsave('~/thesis_figs_and_tables/alp/clf_alp1_clf_control_rescue_donut.pdf', width=10, height=10)

clf_down <- clf %>% filter(FC_KO < 1)
clf_alp2_clf_down <- clf_alp2_clf %>% filter(feature %in% clf_down$feature) %>% dplyr::select(FC_KO, FDR_KO) %>% filter(FC_KO >1)
df <- data.frame(
  Status = c("H3K27me3 Rescued", "H3K27me3 Not Rescued"),
  Count = c(nrow(clf_alp2_clf_down), nrow(clf_down)-nrow(clf_alp2_clf_down))
)
df$ratio <- round((df$Count/sum(df$Count))*100)
p <- donut(df)
ggsave('~/thesis_figs_and_tables/alp/clf_alp2_clf_control_rescue_donut.pdf', width=10, height=10)

##ref6 stats
get_differential_peak_stats('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_midpeak_fixed.csv','/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/')

plot

##Looking at the clf alp2 exclusive mutants from 1.5x FC overlaps
binding_sites <- read.csv('/Volumes/sesame//ALP_Omics/ext_datasets//Shu2019_CLF_SWN_targets/shu_binding_sites.csv')
clf_up <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv') %>% filter(Type.of.K27.reduction.pattern=='III') %>%
  dplyr::select(clf_K27.increase) %>% filter(grepl('AT',clf_K27.increase))
clf_down <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv') %>% filter(Type.of.K27.reduction.pattern=='III') %>%
  dplyr::select(clf_K27.decrease)%>% filter(grepl('AT',clf_K27.decrease))
swn_up <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv') %>% filter(Type.of.K27.reduction.pattern=='III') %>%
  dplyr::select(swn_K27.increase) %>% filter(grepl('AT',swn_K27.increase))
swn_down <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv') %>% filter(Type.of.K27.reduction.pattern=='III') %>%
  dplyr::select(swn_K27.decrease)%>% filter(grepl('AT',swn_K27.decrease))

bound_list <- list(CLF=binding_sites$CLF.bound.genes, SWN=binding_sites$SWN.bound.genes)
binding_venn <- ggVennDiagram::ggVennDiagram(bound_list) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='All Bound Sites') +  scale_color_brewer(palette = "Paired")+
    theme(text = element_text(size = 50),plot.title = element_text(size = 16),legend.title = element_text(size = 12))+
    guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+theme(plot.margin = margin(0, 0, 0, 0))
ggsave('~/thesis_figs_and_tables/alp/SWN_CLF_binding_sites_venn.pdf',binding_venn,height=7, width=7)
K27me3 <- list('clf H3K27me3 Up'=clf_up$clf_K27.increase, 'swn H3K27me3 Up'=swn_up$swn_K27.increase,
'clf H3K27me3 Down'=clf_down$clf_K27.decrease, 'swn H3K27me3 Down'=swn_down$swn_K27.decrease)
me3_venn <- ggVennDiagram::ggVennDiagram(K27me3) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='All H3K27me3 2x Up/Down Peaks') +  scale_color_brewer(palette = "Paired")+
    theme(text = element_text(size = 50),plot.title = element_text(size = 16),legend.title = element_text(size = 12))+
    guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+theme(plot.margin = margin(0, 0, 0, 0))
ggsave('~/thesis_figs_and_tables/alp/SWN_CLF_H3K27me3_sites_venn.pdf',me3_venn,height=7, width=7)
K27me3 <- list('clf H3K27me3 changes'=append(clf_down$clf_K27.decrease,clf_up$clf_K27.increase), 
'swn H3K27me3 changes'=append(swn_up$swn_K27.increase,swn_down$swn_K27.decrease), CLF=binding_sites$CLF.bound.genes, SWN=binding_sites$SWN.bound.genes)
me3_venn <- ggVennDiagram::ggVennDiagram(K27me3) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='All H3K27me3 2x Up/Down Peaks') +  scale_color_brewer(palette = "Paired")+
    theme(text = element_text(size = 50),plot.title = element_text(size = 16),legend.title = element_text(size = 12))+
    guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+theme(plot.margin = margin(0, 0, 0, 0))
ggsave('~/thesis_figs_and_tables/alp/SWN_CLF_H3K27me3_sites_and_H3K27me3_changes_venn.pdf',me3_venn,height=10, width=10)




sicer_swn_up <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_swn_v_col/SRR6453483_sorted-and-SRR6453477_sorted-W200-G600-summary_annotated.csv') %>%
filter(FDR_SWN_vs_Col < 0.05 | FDR_Col_vs_SWN < 0.05) %>% filter(Fc_SWN_vs_Col > 1.2)
sicer_swn_down <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_swn_v_col/SRR6453483_sorted-and-SRR6453477_sorted-W200-G600-summary_annotated.csv') %>%
filter(FDR_SWN_vs_Col < 0.05 | FDR_Col_vs_SWN < 0.05) %>% filter(Fc_SWN_vs_Col <0.8)
swn_list <- list(bound=binding_sites$SWN.unique.bound.genes,swn_up=swn_up$swn_K27.increase, swn_down=swn_down$swn_K27.decrease, sicer_swn_up=sicer_swn_up$feature, sicer_swn_down=sicer_swn_down$feature)


get_unique_comparisons_for_venn <- function(data) {
  unique_values <- unique(names(data))
  num_values <- length(unique_values)
  
  comparisons <- list()
  
  for (i in 1:(num_values - 1)) {
    for (j in (i + 1):num_values) {
      comparisons[[length(comparisons) + 1]] <- c(unique_values[i], unique_values[j])
    }
  }
  
  return(comparisons)
}

empty_df <- data.frame('Set 1' = character(), 'Set 2' = character(), 'pvalue' = numeric(),check.names = FALSE)
comp <- get_unique_comparisons(swn_list)

for(i in 1:length(comp)){
  
  current_comp <- comp[i]
  current_comp_1 <- current_comp[[1]][1]
  current_comp_2 <- current_comp[[1]][2]
  A <- swn_list[names(swn_list)==current_comp[[1]][1]]
  A_filtered <- A[[1]][A[[1]] != ""]
  A <- A_filtered[grep("^AT", A_filtered)]
  B <- swn_list[names(swn_list)==current_comp[[1]][2]]
  B <- B[B!= "-"]
  B_filtered <- B[[1]][B[[1]] != ""]
  B <- B_filtered[grep("^AT", B_filtered)]
  K <- length(intersect(A,B))
  A <- length(setdiff(A,K))
  B <- length(setdiff(B,K))
  print(enrich_pvalue(num_genes_A_specific = A, num_genes_B_specific = B, overlap = K))
  p <- enrich_pvalue(num_genes_A_specific = A, num_genes_B_specific = B, overlap = K)
  outlist <- data.frame('Set 1'= current_comp_1,'Set 2'= current_comp_2, 'pvalue'=p,check.names = FALSE)
  empty_df <- rbind(empty_df, outlist)
}

tab <- ggtexttable(empty_df, rows = NULL, theme = ttheme("blank")) 
tab %>% tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 2, linetype = 2) %>%
tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
 tab_add_hline(at.row = c(nrow(tab)), row.side = "bottom", linewidth = 3, linetype = 1)

overlaps_and_venn_output(swn_list, '~/Desktop/swn_list', plot_title = 'SWN bound sites and swn4 H3K27me3 affected sites',
label_size = 3.5, text_size = 6)
formattable::formattable(empty_df)

l <- find_all_intersections(swn_list)
ggVennDiagram::ggVennDiagram(swn_list) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='All H3K27me3 2x Up/Down Peaks') +  scale_color_brewer(palette = "Paired")+
    theme(text = element_text(size = 50),plot.title = element_text(size = 16),legend.title = element_text(size = 12))+
    guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+theme(plot.margin = margin(0, 0, 0, 0))


swn_all <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_all_binding_and_down/all_swn_bound_and_down.csv')
clf_swn_all <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(feature %in% swn_all$tair) %>%
  mutate(genotype=paste0('clf28 ( n = ', nrow(clf_swn_all),')'))
clf_alp1_swn_all <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_midpeak.csv') %>% filter(feature %in% swn_all$tair) %>%
  mutate(genotype=paste0('clf28 alp1 ( n = ', nrow(clf_swn_all),')'))
clf_alp2_swn_all <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv') %>% filter(feature %in% swn_all$tair) %>%
  mutate(genotype=paste0('clf28 alp2 ( n = ', nrow(clf_swn_all),')'))
swn_all_alp_doubles <- rbind(clf_swn_all, clf_alp1_swn_all, clf_alp2_swn_all)
make_violin_plots_K27(swn_all_alp_doubles,'~/thesis_figs_and_tables/alp/clf_double_mutants_swn_all_sites', neg_num_y =5 , no_neg_num_y = 2.5)


swn_spec <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/swn_down_and_bound_genes.csv')
clf_swn_spec <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(feature %in% swn_spec$x) %>%
  mutate(genotype=paste0('clf28 ( n = ', nrow(clf_swn_spec),')'))
clf_alp1_swn_spec <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_midpeak.csv') %>% filter(feature %in% swn_spec$x) %>%
  mutate(genotype=paste0('clf28 alp1 ( n = ', nrow(clf_swn_spec),')'))
clf_alp2_swn_spec <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv') %>% filter(feature %in% swn_spec$x) %>%
  mutate(genotype=paste0('clf28 alp2 ( n = ', nrow(clf_swn_spec),')'))
swn_spec_alp_doubles <- rbind(clf_swn_spec, clf_alp1_swn_spec, clf_alp2_swn_spec)
make_violin_plots_K27(swn_spec_alp_doubles,'~/thesis_figs_and_tables/alp/clf_double_mutants_swn_specific_sites', neg_num_y = 3.5, no_neg_num_y = 1)


swn_spec <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/swn_down_and_bound_genes.csv')
clf_alp1_swn_spec <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv') %>% filter(feature %in% swn_spec$x) %>%
  mutate(genotype=paste0('clf28 alp1 ( n = ', nrow(clf_swn_spec),')'))
clf_alp2_swn_spec <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv') %>% filter(feature %in% swn_spec$x) %>%
  mutate(genotype=paste0('clf28 alp2 ( n = ', nrow(clf_swn_spec),')'))
swn_spec_alp_doubles <- rbind(clf_swn_spec, clf_alp1_swn_spec, clf_alp2_swn_spec)
make_violin_plots_K27(swn_spec_alp_doubles,'~/thesis_figs_and_tables/alp/clf_double_mutants_swn_specific_sites', neg_num_y = 3.5, no_neg_num_y = 1)



###Transcriptomics tables and overlaps
overall_counts <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/summary/reference_Col-0/summary_DE_genes.csv') %>%
dplyr::filter(grepl('reference_Col-0', Comparison)) %>% filter(grepl('Bennett', Comparison)) %>% filter(grepl('noshrink', Comparison)) %>%
mutate(Comparison=gsub('.csv','', basename(Comparison))) %>%
mutate(Comparison=gsub('_',' ',Comparison)) %>% 
filter(!grepl('hdp', Comparison)) %>% filter(!grepl('ref', Comparison)) %>% dplyr::rename(Genotype=Comparison)

f <- formattable(overall_counts, 
            align =c(rep('c', ncol(overall_counts))), list(Genotype = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))

export_formattable(f, '~/thesis_figs_and_tables/alp/rna/BT_counts.png')

clf_rna_up <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28.csv') %>%
filter(log2FoldChange > 0) %>% dplyr::select(gene_id, log2FoldChange,padj) %>% mutate(genotype='clf28')
write.csv(clf_rna_up, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/broad_analysis/clf_up.csv',
 quote=F,row.names=F)
clf_rna_down <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28.csv') %>%
filter(log2FoldChange < 0)%>% dplyr::select(gene_id, log2FoldChange,padj)%>% mutate(genotype='clf28')
write.csv(clf_rna_down, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/broad_analysis/clf_down.csv',
 quote=F,row.names=F)
clf_alp1_rna_up <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28_alp1-1.csv') %>%
filter(log2FoldChange > 0)%>% dplyr::select(gene_id, log2FoldChange,padj)%>% mutate(genotype='clf28 alpl1')
write.csv(clf_alp1_rna_up, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/broad_analysis/clf_alp1_up.csv',
 quote=F,row.names=F)
clf_alp1_rna_down <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28_alp1-1.csv') %>%
filter(log2FoldChange < 0)%>% dplyr::select(gene_id, log2FoldChange,padj)%>% mutate(genotype='clf28 alp1')
write.csv(clf_alp1_rna_down, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/broad_analysis/clf_alp1_down.csv',
 quote=F,row.names=F)
clf_alp2_rna_up <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28_alp1-1.csv') %>%
filter(log2FoldChange > 0)%>% dplyr::select(gene_id, log2FoldChange,padj)%>% mutate(genotype='clf28 alp2')
write.csv(clf_alp2_rna_up, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/broad_analysis/clf_alp2_up.csv',
 quote=F,row.names=F)
clf_alp2_rna_down <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28_alp1-1.csv') %>%
filter(log2FoldChange < 0)%>% dplyr::select(gene_id, log2FoldChange,padj)%>% mutate(genotype='clf28 alp2')
write.csv(clf_alp2_rna_down, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/broad_analysis/clf_alp2_down.csv',
 quote=F,row.names=F)


# all_fpkms <- read.csv('~/thesis_figs_and_tables/alp/rna/FPKMS_2019.csv')
# df_wide <- as.data.frame(pivot_wider(all_fpkms, names_from = exp, values_from = c(mean_FPKM, se_FPKM)) %>% dplyr::rename(feature=gene_id))
# rownames(df_wide) <- df_wide$feature
# df_wide <- df_wide[grepl('mean', names(df_wide))]
# df_wide <- as.matrix(df_wide)
# scaled <- as.data.frame(t(scale(t(df_wide))))
# write.csv(scaled,'~/thesis_figs_and_tables/alp/rna/scaled_fpkm_mean_values.csv', quote=F)


# all_me3 <- list.files('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params', pattern='me3_annotated_midpeak', full.names = T)
# all_me3 <- list.files('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_macs3', pattern='annotated.csv', full.names = T)
# dataframe <-  data.frame('feature'=character(), 'Score'=numeric(), exp=character())
# for(i in all_me3){
#   title <- (gsub('_annotated_midpeak.csv','',basename(i)))
#   title <- (gsub('_annotated.csv','',basename(i)))
#   df <- read.csv(i)  %>% dplyr::select(feature, Score)%>% mutate(exp=title)
#   print(head(df))
#   dataframe <- rbind(dataframe, df)
# }

# df <- dataframe %>%
#   group_by(feature, exp) %>%
#   summarise(Score = mean(Score))
# me3_df_wide <- as.data.frame(pivot_wider(df, names_from = exp, values_from = c(Score)))
# rownames(me3_df_wide) <- me3_df_wide$feature
# me3_df_wide <- me3_df_wide %>% dplyr::select(!feature)
# me3_df_wide <- as.matrix(me3_df_wide)
# me3_scaled <- as.data.frame(t(scale(t(me3_df_wide)))) %>% na.omit()
# write.csv(me3_scaled,'~/thesis_figs_and_tables/alp/rna/scaled_me3_mean_values.csv')

me3_scaled <- read.csv('~/thesis_figs_and_tables/alp/rna/scaled_me3_mean_values.csv') %>% tibble::column_to_rownames('X') %>% 
dplyr::select(clf, clf.alp1, clf.alp2, Col.0)
rna_scaled <- read.csv('~/thesis_figs_and_tables/alp/rna/scaled_fpkm_mean_values.csv') %>% tibble::column_to_rownames('X') %>% 
dplyr::select(clf, clf.alp1, clf.alp2, Col.0)

swn_down_and_bound <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/swn_down_and_bound_genes.csv')
clf_up_genes <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO >2)
clf_down_genes <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO < 0.5)
clf_any_diff <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv')
me3_scaled_clf_down <- me3_scaled %>% filter(rownames(me3_scaled) %in% clf_down_genes$feature)
scaled_clf_down <- rna_scaled %>% filter(rownames(rna_scaled) %in% rownames(me3_scaled_clf_down))
me3_scaled_swn_down_and_bound <- me3_scaled %>% filter(rownames(me3_scaled) %in% swn_down_and_bound$x)
scaled_swn_down_and_bound <- rna_scaled %>% filter(rownames(rna_scaled) %in% rownames(me3_scaled_swn_down_and_bound))
me3_scaled_clf_up <- me3_scaled %>% filter(rownames(me3_scaled) %in% clf_up_genes$feature)
scaled_clf_up <- rna_scaled %>% filter(rownames(rna_scaled) %in% rownames(me3_scaled_clf_up))
me3_scaled_any_diff <- me3_scaled %>% filter(rownames(me3_scaled) %in% clf_any_diff$feature)
scaled_any_diff <- rna_scaled %>% filter(rownames(rna_scaled) %in% rownames(me3_scaled_any_diff))

p_me3 <- pheatmap::pheatmap(me3_scaled_swn_down_and_bound,treeheight_row = 0, color=colorRampPalette(c("#A6CEE3", "white", "#FF7F00"))(15),
              , fontsize_row = 7,fontsize_col = 10,
              cellheight = 1, cluster_cols=F,cluster_rows=F, angle_col = 90, show_rownames=FALSE,
              cellwidth=10)
p_rna <- pheatmap::pheatmap(scaled_swn_down_and_bound,treeheight_row = 0, color=colorRampPalette(c("#A6CEE3", "white", "#FF7F00"))(15),
              , fontsize_row = 7,fontsize_col = 10,
              cellheight = 1, cluster_cols=F,cluster_rows=F, angle_col = 90, show_rownames=FALSE,
              cellwidth=10)

ggsave('~/Desktop/test_plot.pdf', p_me3, height=100, width=50, unit='mm')
ggsave('~/Desktop/tes_plot.pdf', p_rna, height=100, width=50, unit='mm')






clf_alp1_clf_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv') %>% 
  left_join(df_wide) %>% dplyr::select(feature, FC_KO, mean_FPKM_clf28,"mean_FPKM_clf28 alp2", "mean_FPKM_clf28 alp1", "mean_FPKM_clf28 alp2")
dds <- DESeq(dds)
vsd <- assay(vst(dds))
buds_z <- as.data.frame(t(scale(t(vsd))))

full <- rbind(clf_rna_up, clf_rna_down, clf_alp1_rna_up, clf_alp1_rna_down, clf_alp2_rna_up, clf_alp2_rna_down)
write.csv(full, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/broad_analysis/all_clf_doubles.csv', quote=F, row.names=F)
###Final tables
clf28_stats <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/2019_col_control_peak_stats_clf28_me3_diff_annotated_midpeak.csv')
clf29_stats <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_clf_v_col_merged_bams/peak_stats_clf_merged_me3-and-col_merged_me3-W200-G600-summary_unedited')
swn_stats <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_swn_v_col_merged_bams/peak_stats_swn_merged_me3-and-col_merged_me3-W200-G600-summary_unedited')
all_swn_clf_stats <- rbind(clf28_stats, clf29_stats, swn_stats)
t <- kbl(all_swn_clf_stats, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down","add_linespace=-1.5mm")) %>% 
  row_spec(0, bold = FALSE) %>% 
  column_spec(1, bold = FALSE) %>% 
  kable_styling() %>%
  column_spec(1:ncol(all_swn_clf_stats), border_left = FALSE, border_right = FALSE)
cat(t)


rna_table <- read.csv('~/thesis_figs_and_tables/alp/rna/DEGs_table.csv')
rna_table <- data.frame(lapply(rna_table, function(x) as.character(gsub("_", " ", x))))
t <- kbl(rna_table, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = FALSE) %>% 
  column_spec(1, bold = TRUE) %>% 
  column_spec(1, italic = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:ncol(rna_table), border_left = FALSE, border_right = FALSE)

all_GO <- read.table('~/Salba_RNA/genelists/tair_pantherGO_terms.csv', sep='\t', header = T) %>% dplyr::rename(tair=V2)
true_rescues <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/all_rescues/gene_list.csv') %>% 
dplyr::rename(tair=ensembl_gene_id, 'Gene Name'=external_gene_name) %>% arrange(tair)


num_rows <- nrow(true_rescues)

# Split the dataframe into three parts
df1 <- true_rescues[1:(num_rows/3), ]
df2 <- true_rescues[(num_rows/3 + 1):(2*num_rows/3), ]
df3 <- true_rescues[(2*num_rows/3 + 1):num_rows, ]
df1 <- cbind(df1, df2, df3)
t <- kbl(df1, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = FALSE) %>% 
  column_spec(c(1,3,5), bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:ncol(df1), border_left = FALSE, border_right = FALSE)



####Final venns
##BT-clf Shu et al
shu_peaks_swn <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv') %>% dplyr::select(swn_K27.decrease) %>% filter(grepl('AT',swn_K27.decrease ))
shu_peaks_clf <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv') %>% dplyr::select(clf_K27.decrease)%>% filter(grepl('AT',clf_K27.decrease ))
genelist= list('swn-4'= shu_peaks_swn$swn_K27.decrease, 'clf-29'=shu_peaks_clf$clf_K27.decrease)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/swn4_clf29_sicer_comparison_2x_down',
stats='no',plot_title='Overlap of Differential H3K27me3 relative to Col-0 in the two datasets', text_size = 7.5)


bt <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv')
shu <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_clf_v_col_merged_bams/clf_merged_me3-and-col_merged_me3-W200-G600-summary_annotated.csv')
genelist <- list( 'Shu_clf-29'=shu$feature,'Bennet_clf-28'=bt$feature)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/clf28_clf29_sicer_comparison_2x',
stats='no',plot_title='Overlap of Differential H3K27me3 relative to Col-0 in the two datasets', text_size = 7.5)

##SHU CLF-SWN 
clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_clf_v_col_merged_bams/clf_merged_me3-and-col_merged_me3-W200-G600-summary_annotated.csv') %>%
filter(Fc_CLF_vs_Col > 1.5 |  Fc_CLF_vs_Col < 0.75)
swn <-read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_swn_v_col_merged_bams/swn_merged_me3-and-col_merged_me3-W200-G600-summary_annotated.csv') %>%
filter(Fc_SWN_vs_Col > 1.5 |  Fc_SWN_vs_Col < 0.75)
genelist <- list( 'Shu_clf-29'=clf$feature,'Shu_swn-4'=swn$feature)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/clf29_swn_sicer_comparison_1.2x',
stats='no',plot_title='Overlap of Differential H3K27me3 relative to Col-0 in the two datasets', text_size = 7.5)

##SHU CLF-SWN Bound
clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_clf_v_col_merged_bams/clf_merged_me3-and-col_merged_me3-W200-G600-summary_annotated.csv') %>%
filter(Fc_CLF_vs_Col > 1.5 |  Fc_CLF_vs_Col < 0.75)
swn <-read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_swn_v_col_merged_bams/swn_merged_me3-and-col_merged_me3-W200-G600-summary_annotated.csv') %>%
filter(Fc_SWN_vs_Col > 1.5 |  Fc_SWN_vs_Col < 0.75)
clf_bound <- readxl::read_xlsx('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/PLD3-3-e00100-s004_Bound_Genes.xlsx') %>% dplyr::select(`CLF-bound genes`)
swn_bound <- readxl::read_xlsx('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/PLD3-3-e00100-s004_Bound_Genes.xlsx') %>% dplyr::select(`SWN-bound genes`)
genelist <- list( 'clf-29 H3K27me3'=clf$feature,'swn-4 H3K27me3'=swn$feature, 'CLF-Bound'=clf_bound$`CLF-bound genes`,'SWN-Bound'=swn_bound$`SWN-bound genes`)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/clf29_swn_sicer_comparison_1.5x_and_bound',
,plot_title='Overlap of Differential H3K27me3 relative to Col-0 in the two datasets', text_size = 5, label_size=5)
go_dot

##clf alp doubles
clfalp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1.2 | FC_KO < 0.8)
clfalp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1.2 | FC_KO < 0.8)
genelist <- list('clf alp1'= clfalp1$feature,'clf alp2'=clfalp2$feature)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/clf_alp_double_mutants_col_control_overlap',
plot_title='Overlap of Differential H3K27me3 relative to Col-0 in the two datasets', stats='no', text_size=10)

clfalp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.2 | FC_KO < 0.8)
clfalp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.2 | FC_KO < 0.8)
genelist <- list('clf alp1'= clfalp1$feature,'clf alp2'=clfalp2$feature)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/clf_alp_double_mutants_clf_control_overlap',
plot_title='Overlap of Differential H3K27me3 relative to clf-28 in the two datasets', stats='no', text_size=10)


#BUD DEG plot
clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv')
clfalp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_midpeak.csv')
clfalp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv')
samps <- c('clf', 'clf alp1','clf alp2')
condition <- c(rep("H3K27me3 Up" , 3),rep("H3K27me3 Down",3))
value <- c(nrow(clf[clf$FC_KO>1,]),nrow(clfalp1[clfalp1$FC_KO>1,]),nrow(clfalp2[clfalp2$FC_KO>1,]),
nrow(clf[clf$FC_KO<1,]),nrow(clfalp1[clfalp1$FC_KO<1,]),nrow(clfalp2[clfalp2$FC_KO<1,]))
data <- data.frame(samps,condition,value)
p_all <- ggplot(data, aes(fill=condition, y=value, x=samps,label = value)) + 
      theme +
      geom_bar(position="stack", stat="identity") +
      geom_text(size = 5, position = position_stack(vjust = 0.5),colour = "white")+
      scale_fill_manual(values=c("H3K27me3 Down" = "#E1BE6A", "H3K27me3 Up" = "#40B0A6"))+
      scale_x_discrete(limits = c("clf", "clf alp1", "clf alp2"))+
      ggtitle("Comparison of H3K27me3 relative to Col-0 - All peaks") +
      theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 20),
      axis.ticks.x = element_blank(),  # Remove ticks
      panel.grid.major.x = element_blank(),  # Remove gridlines
      panel.grid.minor.x = element_blank(),  # Remove gridlines
      legend.title = element_text(size = 15),  # Change legend title and size
      legend.text = element_text(size = 15))+
      guides(fill=guide_legend(title="DEG Status"))+
      xlab("Genotype")+
      ylab("Peak count")

p_1.5 <- ggplot(data, aes(fill=condition, y=value, x=samps,label = value)) + 
      theme +
      geom_bar(position="stack", stat="identity") +
      geom_text(size = 5, position = position_stack(vjust = 0.5),colour = "white")+
      scale_fill_manual(values=c("H3K27me3 Down" = "#E1BE6A", "H3K27me3 Up" = "#40B0A6"))+
      scale_x_discrete(limits = c("clf", "clf alp1", "clf alp2"))+
      ggtitle("Comparison of H3K27me3 relative to Col-0 - 1.5x FC") +
      theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 20),
      axis.ticks.x = element_blank(),  # Remove ticks
      panel.grid.major.x = element_blank(),  # Remove gridlines
      panel.grid.minor.x = element_blank(),  # Remove gridlines
      legend.title = element_text(size = 15),  # Change legend title and size
      legend.text = element_text(size = 15))+
      guides(fill=guide_legend(title="DEG Status"))+
      xlab("Genotype")+
      ylab("Peak count")



p_2 <-  ggplot(data, aes(fill=condition, y=value, x=samps,label = value)) + 
      theme +
      geom_bar(position="stack", stat="identity") +
      geom_text(size = 5, position = position_stack(vjust = 0.5),colour = "white")+
      scale_fill_manual(values=c("H3K27me3 Down" = "#E1BE6A", "H3K27me3 Up" = "#40B0A6"))+
      scale_x_discrete(limits = c("clf", "clf alp1", "clf alp2"))+
      ggtitle("Comparison of H3K27me3 relative to Col-0 - 2x FC") +
      theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 20),
      axis.ticks.x = element_blank(),  # Remove ticks
      panel.grid.major.x = element_blank(),  # Remove gridlines
      panel.grid.minor.x = element_blank(),  # Remove gridlines
      legend.title = element_text(size = 15),  # Change legend title and size
      legend.text = element_text(size = 15))+
      guides(fill=guide_legend(title="DEG Status"))+
      xlab("Genotype")+
      ylab("Peak count")




ggsave('~/thesis_figs_and_tables/alp/fold_change_barplots/clf_alp_double_mutants_all_FC_barplots_combined.pdf', d, height=10, width=10)

d <- ggpubr::ggarrange(p_all,p_1.5,p_2, nrow=1, ncol=3,common.legend = TRUE, legend="bottom")
##clf alp RNA
clf_chip <-read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO >1.5 | FC_KO < 0.75)
clf_rna <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28.csv')
genelist <- list('clf ChIP-seq'= clf_chip$feature,'clf RNA-seq'=clf_rna$gene_id)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/clf_chip_rna_no_fc',
plot_title='Overlap of Differential H3K27me3 relative to clf-28 in the two datasets', stats='no', text_size=8)


clf_alp1_chip <-read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO >1.5 | FC_KO < 0.75)
clf_alp1_rna <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28_alp1-1.csv')
genelist <- list('clf alp1 ChIP-seq'= clf_alp1_chip$feature,'clf alp1 RNA-seq'=clf_alp1_rna$gene_id)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/clf_alp1_chip_rna_no_fc',
plot_title='Overlap of Differential H3K27me3 relative to clf-28 in the two datasets', stats='no', text_size=8)

clf_alp2_chip <-read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO >1.5 | FC_KO < 0.75)
clf_alp2_rna <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28_alp2-1.csv')
genelist <- list('clf alp2 ChIP-seq'= clf_alp1_chip$feature,'clf alp2 RNA-seq'=clf_alp1_rna$gene_id)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/clf_alp2_chip_rna_no_fc',
plot_title='Overlap of Differential H3K27me3 relative to clf-28 in the two datasets', stats='no', text_size=8)

clf_down_sig <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO < 0.75) 
clf_alp2_corresonding <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv') %>% filter(feature %in% clf_down_sig$feature)
df <- data.frame('Status'= c('H3K27me3 Up', 'H3K27me3 Down'),
'Count'=c(nrow(clf_alp2_corresonding[clf_alp2_corresonding$FC_KO >1,]),nrow(clf_alp2_corresonding[clf_alp2_corresonding$FC_KO <1,]))) %>%
mutate(ratio=round((Count/cumsum(Count)[2])*100),digits = 1)
clf_alp2 <- donut(df=df)
ggsave('~/thesis_figs_and_tables/alp/venn_final/clf_alp2_donut_rescue.pdf',clf_alp2)
clf_alp1_corresonding <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv') %>% filter(feature %in% clf_down_sig$feature)
df <- data.frame('Status'= c('H3K27me3 Up', 'H3K27me3 Down'),
'Count'=c(nrow(clf_alp1_corresonding[clf_alp1_corresonding$FC_KO >1,]),nrow(clf_alp1_corresonding[clf_alp1_corresonding$FC_KO <1,]))) %>%
mutate(ratio=round((Count/cumsum(Count)[2])*100),digits = 1)
clf_alp1 <- donut(df=df)
ggsave('~/thesis_figs_and_tables/alp/venn_final/clf_alp1_donut_rescue.pdf',clf_alp1, height=10, width=10)


##ref6 col-0 checks
lu <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/lu_et_al_sicer/ref_me3_sorted-and-col_me3_sorted-W200-G600-summary_annotated.csv')
cui <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/cui_et_al_sicer/ref6_me3_sorted-and-col_me3_sorted-W200-G600-summary_annnotated.csv')
cv <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/sicer_ref6_v_UDIcol_merged_bams/ref6-and-col2022_me3-W200-G600-summary_annnotated.csv.csv')%>%
filter(feature %in% lu$feature) %>% filter(feature %in% cui$feature)
all_genelist <- list('Lu et al 2011'=lu$feature,'Cui et al 2016'=cui$feature, '2021 ref6'=cv$feature)
overlaps_and_venn_output(genes=all_genelist,out_file_prefix='~/thesis_figs_and_tables/alp/ref6/ref6_col_control_2021_peaks_overlap_with_ext_annotated_2xFC',
plot_title='Overlap of differential H3K27me3 sites in ref6 relative to Col-0', text_size=5, label_size=5)

ref6_alp2_ref_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_midpeak.csv') %>% filter(FC_KO >1.5)
clf_alp2_clf_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv') %>% filter(FC_KO >1.5)
genelist <- list('ref6 alp2 - ref6 control'= ref6_alp2_ref_control$feature, 'clf alp2 - clf control'=clf_alp2_clf_control$feature)
overlaps_and_venn_output(genes=genelist, out_file_prefix='~/thesis_figs_and_tables/alp/venn_final/ref6_alp2_ref_control_clf_alp2_clf_control_overlap',
plot_title='Differential H3K27me3 sites in ref6 alp2 relative to ref6 and clf alp2 relative to clf', text_size = 5)
##### GO plots
shu_CLF <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Shu_sicer/sicer_clf_v_col_merged_bams/clf_merged_me3-and-col_merged_me3-W200-G600-summary_annotated.csv') %>%
filter(Fc_CLF_vs_Col > 1.5 |  Fc_CLF_vs_Col < 0.75)
panther_go_maker(shu_CLF$feature,'~/thesis_figs_and_tables/alp/GO_final/clf29_1.5x_panther_GO')
go_dotplotter('~/thesis_figs_and_tables/alp/GO_final/clf29_1.5x_panther_GO')
bt_CLF <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>%
filter(FC_KO > 1.75 |  FC_KO < 0.65)
panther_go_maker(bt_CLF$feature,'~/thesis_figs_and_tables/alp/GO_final/clf28_1.5x_panther_GO')
go_dotplotter('~/thesis_figs_and_tables/alp/GO_final/clf28_1.5x_panther_GO')


write.csv(cv, '/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/sicer_ref6_v_UDIcol_merged_bams/ref6-and-col2022_me3-W200-G600-summary_annnotated_ext_overlap', quote=F, row.names=F)


get_differential_peak_stats_sicer('/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/sicer_ref6_v_UDIcol_merged_bams/ref6-and-col2022_me3-W200-G600-summary_annnotated_ext_overlap',
'/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/sicer_ref6_v_UDIcol_merged_bams/')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/sicer_ref6_v_UDIcol_merged_bams/peak_stats_ref6-and-col2022_me3-W200-G600-summary_annnotated_ext_overlap',
outname='~/thesis_figs_and_tables/alp/fold_change_barplots/ref6_col_control_overlap_with_external_fold_change.png')

get_differential_peak_stats_sicer('/Volumes/sesame/ALP_Omics/ChIP/ref6_me1_sicer/sicer_ref6_alp2_me1_ref6_control/ref6_alp2_me1-and-ref6_me1-W200-G600-summary_annotated.csv',
'/Volumes/sesame/ALP_Omics/ChIP/ref6_me1_sicer/sicer_ref6_alp2_me1_ref6_control/')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/ref6_me1_sicer/sicer_ref6_alp2_me1_ref6_control/peak_stats_ref6_alp2_me1-and-ref6_me1-W200-G600-summary_annotated.csv', 
outname='~/thesis_figs_and_tables/alp/fold_change_barplots/ref6_alp2_me1_ref6_control_fold_changes.png', control='ref6')

panther_go_maker(cv$feature, '/Volumes/sesame/ALP_Omics/ChIP/ref6_sicer/sicer_ref6_v_UDIcol_merged_bams/panther_GO_overlaps_no_fc')
get_differential_peak_stats('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_midpeak_fixed.csv',
'/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/')
plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/peak_stats_ref6_alp2_me3_diff_ref6_control_annotated_midpeak_fixed.csv',
'~/thesis_figs_and_tables/alp/fold_change_barplots/ref6_alp2_ref6_control_fold_changes.png', control='ref6')
#####Final plot figures
peaks_of_interest <- read.csv('~/thesis_figs_and_tables/alp/peak_figures_final/peaks_of_interest.csv')
##CLF SWN - Shu BT comparison
path <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/'
shu_bigwigs <- list('/Volumes/sesame/joerecovery/Project_folder/alp_omics/shu_et_al_2019_bigwigs/bigwigs/clf_me3_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/shu_et_al_2019_bigwigs/bigwigs/swn_me3_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/shu_et_al_2019_bigwigs/bigwigs/col-0_me3_rpkm.bw')
shu_titles <- c('clf', 'swn','Col-0')
plot_peaks(bigwig_list=shu_bigwigs, chrom=3,genename='AT3G23130',start=8242201-2000, end=8243427+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT3G23130_swn_clf_compensation', title_list = shu_titles)
plot_peaks(bigwig_list=shu_bigwigs, chrom=1,genename='AT1G41830',start=15602686-2000, end=15608844+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT1G41830_swn_clf_compensation', title_list = shu_titles)
plot_peaks(bigwig_list=shu_bigwigs, chrom=2,genename='AT2G23510',start=10006992-2000, end=10019093+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT2G23510_swn_clf_compensation', title_list = shu_titles)
plot_peaks(bigwig_list=shu_bigwigs, chrom=5,genename='AT5G41250',start=16500879-2000, end=16503238+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT5G41250_swn_slightly_lower', title_list = shu_titles)
plot_peaks(bigwig_list=shu_bigwigs, chrom=2,genename='AT2G15020',start=6491006-2000, end=6493929+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT2G15020_swn_slightly_lower', title_list = shu_titles)
plot_peaks(bigwig_list=shu_bigwigs, chrom=4,genename='AT4G18960',start=1038000-100, end=10390000+100,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT4G18960_swn_slightly_lower', title_list = shu_titles)
plot_peaks(bigwig_list=shu_bigwigs, chrom=1,genename='AT1G59970',start=22075000-2000, end=22085000+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT1G59970_clf_slightly_lower', title_list = shu_titles)

##alp singles
alp_single_bigwigs <- list('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2020_ChIP_bams/alp1_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2020_ChIP_bams/alp2_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2020_ChIP_bams/col_rpkm.bw')
alp_single_titles <- c('alp1', 'alp2','Col-0')
plot_peaks(bigwig_list=alp_single_bigwigs, chrom=3,genename='AT3G23130',start=8242201-2000, end=8243427+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT3G23130_alp_singles', title_list = alp_single_titles)
plot_peaks(bigwig_list=alp_single_bigwigs, chrom=5,genename='AT5G56420',start=22762000-2000, end=22762000+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT5G56420_alp_singles_outliers', title_list = alp_single_titles)
plot_peaks(bigwig_list=alp_single_bigwigs, chrom=3,genename='AT3G21220',start=7446000-2000, end=7446000+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT3G21220_alp_singles_outliers', title_list = alp_single_titles)
plot_peaks(bigwig_list=alp_single_bigwigs, chrom=1,genename='AT1G63880',start=23712600-2000, end=23715999+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT1G63880_alp_singles', title_list = alp_single_titles)
plot_peaks(bigwig_list=alp_single_bigwigs, chrom=2,genename='AT2G28725',start=12327000-2000, end=12328999+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT2G28725_alp_singles', title_list = alp_single_titles)

##clf alp doubles
clf_alp_bigwigs <- list('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/col_rpkm.bw')
clf_alp_titles <- c('clf', 'clf alp1','clf alp2','Col-0')
plot_peaks(bigwig_list=clf_alp_bigwigs, chrom=4,genename='AT4G18960',start=10380000-2000, end=10390000+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT4G18960_clf_alp_rescue', title_list = clf_alp_titles)
plot_peaks(bigwig_list=clf_alp_bigwigs, chrom=5,genename='AT5G10140',start=3172167, end=3180659,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT5G10140_clf_alp_rescue', title_list = clf_alp_titles)
plot_peaks(bigwig_list=clf_alp_bigwigs, chrom=2,genename='AT2G22800',start=9704315-1000, end=9706534,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT2G22800_clf_alp_rescue', title_list = clf_alp_titles)
plot_peaks(bigwig_list=clf_alp_bigwigs, chrom=5,genename='AT5G22570',start=7495455-2000, end=7496990+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT5G22570_clf_alp_rescue', title_list = clf_alp_titles)
plot_peaks(bigwig_list=clf_alp_bigwigs, chrom=1,genename='AT1G08860',start=611794-2000, end=612923+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT1G08860_clf_alp_rescue', title_list = clf_alp_titles)
plot_peaks(bigwig_list=clf_alp_bigwigs, chrom=1,genename='AT1G01010',start=3518-2000, end=6012+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT1G01010_clf_alp', title_list = clf_alp_titles)


ref_bigwigs <- list('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_alp2_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/alp2_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/col_rpkm.bw')
ref_titles <- c('ref6', 'ref6 alp2','alp2','Col-0')
plot_peaks(bigwig_list=ref_bigwigs, chrom=4,genename='AT4G18960',start=10380000-2000, end=10390000+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT4G18960_ref6_top_peaks', title_list = ref_titles)
plot_peaks(bigwig_list=ref_bigwigs, chrom=3,genename='AT3G45750',start=16799000-2000, end=16800399+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT3G45750_ref6_top_peaks', title_list = ref_titles)
plot_peaks(bigwig_list=ref_bigwigs, chrom=4,genename='AT4G16880',start=9497600-2000, end=9499199+2000,
 out='~/thesis_figs_and_tables/alp/peak_figures_final/AT4G16880_ref6_top_peaks', title_list = ref_titles)


##RNA_FPKM plots
all_fpkms <- read.csv('~/thesis_figs_and_tables/alp/rna/FPKMS_2019.csv') %>% dplyr::rename(genotype=exp)
plot_fold_changes(all_fpkms,'AT5G10140',outprefix = '~/thesis_figs_and_tables/alp/fpkm/AT5G10140')
