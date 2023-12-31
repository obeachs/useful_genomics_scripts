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
library(hypergea)
library(gmp)
??export_formattable
source('/Volumes/sesame/joerecovery/scripts/formattable_functions_and_tables.R')
theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)
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
clfalp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 2)
clfalp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 2)
overlaps_and_venn_output(genes=genelist, out_file_prefix = '~/thesis_figs_and_tables/alp/col_control_clf_doubles', plot_title = 'clf28 and clf28 double mutants - col-0 reference')
clfalp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.5| FC_KO < 0.75)
clfalp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.5 | FC_KO < 0.75)
overlaps_and_venn_output(genes=genelist, out_file_prefix = '~/thesis_figs_and_tables/alp/clf_control_clf_doubles', plot_title = 'clf28 double mutants - clf28 reference')

clf_alp1_clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.5)
clf_alp2_clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv')%>% filter(FC_KO > 1.5)
alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1)
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_midpeak.csv')%>% filter(FC_KO > 1)









swn_down_and_bound <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/swn_overlap/SWN_specific_bound_and_down/swn_down_and_bound_genes.csv')

genelist <- list(clf_alp1=clfalp1$feature, clf_alp2=clfalp2$feature,swn=swn_down_and_bound$x)
overlaps_and_venn_output(genes=genelist, out_file_prefix = '~/thesis_figs_and_tables/alp/clf_alp_doubles_2x_UP_FC_swn_down_and_bound_overlaps', plot_title = 'clf alp double mutants and swn', text_size = 10)


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
  scale_fill_manual(values = c("H3K27me3 Down" = "#f96955", "H3K27me3 Up" = "#aff581")) +
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

print(p)



ggsave('~/Desktop/test_plot.pdf',p, width=10, height=10,dpi = 45500)

d <- ggviolin(clf_volcano_all, x = "genotype", y = "FC_KO", palette = "Paired", fill = "genotype") + 
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
    breaks = seq(round(min(clf_volcano_all$FC_KO)) - 1, round(max(clf_volcano_all$FC_KO)) + 1, by = 1),
    limits = c(min(clf_volcano_all$FC_KO),max(clf_volcano_all$FC_KO))
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab('H3K27me3 Fold Change Relative to Col-0') +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")
  # stat_compare_means(
  #   comparisons = my_comparisons,
  #   method = "wilcox.test",  # You can use other methods such as "wilcox.test", "anova", etc.
  #   #label = "p.signif",
  #   #size = 4,
  #   #step.increase = 0.2




ggsave('~/thesis_figs_and_tables/alp/clf28_col_control_H3K27me3_fold_changes_volcanoplot_no_neg.pdf',d, width=10, height=10)
ggsave('~/thesis_figs_and_tables/alp/clf28_col_control_H3K27me3_fold_changes_barplot.pdf',p, width=10, height=10)


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
  scale_fill_manual(values = c("H3K27me3 Down" = "#f96955", "H3K27me3 Up" = "#aff581")) +
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
  scale_fill_manual(values = c("H3K27me3 Down" = "#f96955", "H3K27me3 Up" = "#aff581")) +
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

tab <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats//diff_peaks/2019_clf_control_peak_stats_clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv')
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
    scale_fill_manual(values = c("H3K27me3 Down" = "#f96955", "H3K27me3 Up" = "#aff581")) +
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


  ggsave(outname,p, width=10, height=10)
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

plot_peak_stats_barplot('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/alp2_FC_ref_control_table.csvpeak_stats_alp2_me3_diff_ref6_control_annotated_midpeak.', 
'/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/alp2_FC_ref_control_table.csvpeak_stats_alp2_me3_diff_ref6_control_annotated_midpeak.pdf', control='ref6')

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

full <- rbind(clf_rna_up, clf_rna_down, clf_alp1_rna_up, clf_alp1_rna_down, clf_alp2_rna_up, clf_alp2_rna_down)
write.csv(full, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/broad_analysis/all_clf_doubles.csv', quote=F, row.names=F)
