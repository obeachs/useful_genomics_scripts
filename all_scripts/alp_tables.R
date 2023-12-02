library(ggplot2)
library(ggplot2)
library(dplyr)
library(ggVennDiagram)
library(gplots)
library(VennDiagram)
library(RVenn)
library(stringr)
library(ggpubr)
library(ggrepel)
library(tidyverse)
source('/Volumes/sesame/joerecovery/scripts/formattable_functions_and_tables.R')



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


diff <- list.files('/Volumes/sesame/ALP_Omics//ChIP/validations//stats//peak_number_stats/diff_peaks/',patter='m',full.names = T)
diff_test <- read.csv(diff[1])
diff_df <- matrix(ncol = length(names(diff_test)), nrow = 0)
colnames(diff_df) <- names(diff_test)
for(i in diff){
    df <- read.csv(i)
    rownames(df) <- basename(i)
    diff_df <- rbind(diff_df, df)
}


write.csv(diff_df,'/Volumes/sesame/ALP_Omics/ChIP/validations/stats/peak_number_stats/diff_peaks/all_stats.csv', quote=F, row.names=F)

diff_df_col <- diff_df %>% filter(grepl('2020', rownames(diff_df)) | grepl('col_control', rownames(diff_df))) %>%
rownames_to_column('Genotype') %>% mutate(Genotype=gsub("^.*peak_stats_(.*?)\\..*$", "\\1", Genotype)) %>% 
mutate(Genotype=gsub('_diff_',' ', Genotype)) %>% mutate(Genotype=gsub('_option','', Genotype)) %>% mutate(Genotype=gsub('_fixed',' ', Genotype))%>%
mutate(Genotype=gsub('_',' ', Genotype)) %>% mutate(Genotype=gsub('annotated',' ', Genotype)) %>% mutate(Genotype=gsub('annotations',' ', Genotype)) %>%
mutate(Genotype=gsub('both','(both)', Genotype))  %>% mutate(Genotype=gsub('midpeak','(nearest)', Genotype)) %>%
mutate("Histone Mark"=ifelse(grepl('me3', Genotype),'H3K27me3','H3K27me1')) %>% mutate(Genotype=gsub('me3','', Genotype)) %>%
mutate(Genotype=gsub('me1','', Genotype)) %>% filter(!grepl('transposons', Genotype))
names(diff_df_col) <- gsub('_',' ',names(diff_df_col))
names(diff_df_col) <- str_to_title(names(diff_df_col))

make_nice_table(df=diff_df_col, outname='~/thesis_figs_and_tables/alp/Differential_peaks_table.png')
f <- formattable(diff_df_col, 
            align =c(rep('c', ncol(diff_df_col))), list(Genotype = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))
export_formattable(f,'~/thesis_figs_and_tables/alp/Differential_peaks_table.png')

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




shu <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv', skip=1) %>% filter(Type.of.K27.reduction.pattern=='III')
bt <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(FC_KO >2 | FC_KO < 0.5)
clf_up <-length(shu$clf_K27.increase[grepl('AT',shu$clf_K27.increase)])
clf_down <- length(shu$clf_K27.decrease[grepl('AT',shu$clf_K27.decrease)])
swn_up <- length(shu$swn_K27.increase[grepl('AT',shu$swn_K27.increase)])
swn_down <- length(shu$swn_K27.decrease[grepl('AT',shu$swn_K27.decrease)])
bt_up <- length(bt$feature[bt$FC_KO>1])
bt_down <- length(bt$feature[bt$FC_KO<1])


df <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/H3K27me3_counts.csv')
colnames(df) <- c('Experiment', 'H3K27me3 Up', 'H3K27me3 Down')

f <- formattable(df, 
            align =c(rep('c', ncol(df))), list(Experiemnt = formatter(
              "span", style = ~ style(color = "grey",font.style = "italic"))))
export_formattable(f,'~/thesis_figs_and_tables/alp/shu_bt_overall_H3K27me3.png')


