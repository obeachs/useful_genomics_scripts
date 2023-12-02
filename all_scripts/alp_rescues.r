library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gplots)
library(VennDiagram)
library(RVenn)
library(stringr)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(pals)
library(clusterProfiler)
library(drawProteins)
library(ggpubr)
library(readxl)
library(knitr)
library(org.At.tair.db)
library(rtracklayer)
library('coriell')
library(topGO)
source('~/ChIP_validations.R')
library(biomaRt)
library(httr)
library(XML)
source('/Volumes/sesame/ALP_Omics/ChIP/validations/alp_visualisation_scripts.r')
source('/Volumes/sesame/ALP_Omics/ChIP/scripts/annotate_peaks.R')
source('~/Salba_RNA/scripts/GO_dotplotter.R')


annotate_peaks('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/alp2_me3_diff.txt', chr_col='')
find_all_intersections <- function(list_of_files){
  intersection_total <- gplots::venn(list_of_files,show.plot = FALSE)
  intersection_list <- attributes(intersection_total)$intersections
  intersection_list
}


 mart <- biomaRt::useMart(
    host = "https://plants.ensembl.org", biomart = "plants_mart",
    dataset = "athaliana_eg_gene"
  )




gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   values = rescue_genelist,
                   mart = mart)

print(gene_info$external_gene_name)

rescues <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations//clf_double_mutants/rna_chip_comparison/all_rescues/gene_list.csv') %>%dplyr::select(ensembl_gene_id)
gene_info <- filter(gene_info, ensembl_gene_id %in% rescue_genelist_clf_down)


p <- general_venn(list(res=rescues$ensembl_gene_id, feature= trans_clf_alp2$gene_id))
ggsave('~/Desktop/test_plot.pdf', p)
#write.csv(gene_info$external_gene_name,'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/gene_list_names.csv', quote=F, row.names=F)


trans_clf <-read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28.csv') %>% mutate(exp="clf28")
trans_clfalp1 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28_alp1-1.csv')%>% mutate(exp="clf28_alp1")
trans_clfalp2 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28_alp2-1.csv')%>% mutate(exp="clf28_alp2")

counts <- list.files('/Volumes/sesame/ALP_Omics/RNASeq/counts', full.names = T)
data <- data.frame()
for(i in counts){
  name <- basename(i)
  exp_name <- sub("^(.*)_.*$", "\\1", name)
  print(exp_name)
  df <- read.table(i, header = T) %>% mutate(exp=exp_name)
  print(head(df))
  data <- rbind(df, data)
}

for(i in 1:length(rescues$ensembl_gene_id)){
  plot_fold_changes(data, rescues$ensembl_gene_id[i], paste0'~/thesis_figs_and_tables/alp/')
}

summary_data <- FLC %>%
  group_by(exp) %>%
  summarise(
    mean_FPKM = mean(FPKM),
    sd_FPKM = sd(FPKM),
    se_FPKM = sd_FPKM / sqrt(n())
  )

# Create a bar plot with error bars
ggplot(summary_data, aes(x = exp, y = mean_FPKM, fill = exp)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_FPKM - se_FPKM, ymax = mean_FPKM + se_FPKM), position = position_dodge(0.9), width = 0.25) +
  labs(title = "Average FPKM with Error Bars Grouped by exp",
       x = "exp",
       y = "Average FPKM") +
  theme_minimal()


df <- rbind(trans_clf,trans_clfalp1,trans_clfalp2 )
for(i in rescues$ensembl_gene_id){
  plot_fold_changes(data, i, '~/thesis_figs_and_tables/alp/FPKM')
}

trans_clf <-read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_Col-0/clf-28.csv') %>% dplyr::select(gene_id, log2FoldChange, padj) %>% 
mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05) %>% filter(log2FoldChange > 0)
chip_clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% filter(FC_WT >1.5, FDR_WT<0.05) %>% 
filter(feature %in% trans_clf$gene_id)

clf_alp1_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv') %>% filter(feature %in% chip_clf$feature)
clf_alp2_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv') %>% filter(feature %in% chip_clf$feature)

trans_clf_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp1-1.csv') %>% dplyr::select(gene_id, log2FoldChange, padj) %>% 
mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05) %>% filter(gene_id %in% clf_alp1_control$feature)
trans_clf_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp2-1.csv') %>% dplyr::select(gene_id, log2FoldChange, padj) %>% 
mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05)%>% filter(gene_id %in% clf_alp2_control$feature) %>% filter(log2FoldChange > 0)

trans_clf_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp2-1.csv') %>% filter(gene_id %in% rescues$ensembl_gene_id) %>%
dplyr::select(gene_id, log2FoldChange, padj)
trans_clf_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp1-1.csv') %>% filter(gene_id %in% rescues$ensembl_gene_id) %>%
dplyr::select(gene_id, log2FoldChange, padj)
clf_alp1_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_midpeak_annotated.csv') %>% 
filter(feature %in% rescues$ensembl_gene_id) %>% dplyr::select(feature, FC_KO, FC_WT)
clf_alp2_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv') %>% 
filter(feature %in% rescues$ensembl_gene_id)%>% dplyr::select(feature, FC_KO, FC_WT)
inner_join(trans_clf_alp1, trans_clf_alp2, by='gene_id')

alp2_bt <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/transposon_check/txdb_transposable_elements/alp2_me3_diff_annotations_transposons.csv') %>% filter(FC_KO > 1.5)
alp2_bt_stats <- data.frame(c)

alp2_christos <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/alp2_me3_diff_annotated_midpeak_transposons_fixed.csv') %>% filter(FC_KO > 1.5)
panther_go_maker(intersect(alp2_bt$gene_id, alp2_christos$feature), '~/thesis_figs_and_tables/alp/alp2_2020_me3_v_2021_me3_intesect_FC1.5')
go_dotplotter('~/thesis_figs_and_tables/alp/alp2_2020_me3_v_2021_me3_intesect_FC1.5')
alp2_me3 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/alp2_me3_diff_annotated_midpeak_transposons_fixed.csv') %>% filter(FC_KO > 1.5)
alp2_me1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/alp2_me1_diff_annotated_midpeak_transposons_fixed.csv') %>% filter(FC_KO > 1.5)
p <- general_venn(list(me1=alp2_me1$feature,me3=alp2_me3$feature), title='Overlap of H3K27me1 and H3K27me3 peaks - Up relative to Col-0')
overlaps_and_venn_output(list(me1=alp2_me1$feature,me3=alp2_me3$feature),out_file_prefix = '~/thesis_figs_and_tables/alp/alp2_2021_me1_v_me3_overlaps_FC_UP_overlap.pdf',plot_title='H3K27me1 and H3K27me1 overlap')
ggsave('~/thesis_figs_and_tables/alp/alp2_2021_me1_v_me3_overlaps_FC_up_only.pdf',p)


overlaps_and_venn_output(list(chr2021=alp2_christos$feature,bt2020=alp2_bt$gene_id),out_file_prefix = '~/thesis_figs_and_tables/alp/alp2_2020_me3_v_2021_me3_overlaps_FC_UP_overlap.pdf',plot_title='2020 and 2021 overlap')
panther_go_maker(alp2_bt$feature, '~/thesis_figs_and_tables/alp/alp2_2020_me3_pantherGO')
go_dotplotter('~/thesis_figs_and_tables/alp/alp2_2020_me3_pantherGO')
panther_go_maker(alp2_christos$feature, '~/thesis_figs_and_tables/alp/alp2_2021_me3_pantherGO')
go_dotplotter('~/thesis_figs_and_tables/alp/alp2_2021_me3_pantherGO')

panther_go_maker(intersect(alp2_me3$feature, alp2_me1$feature), '~/thesis_figs_and_tables/alp/alp2_2021_me1_v_me3_overlaps_pantherGO')
go_dotplotter('~/thesis_figs_and_tables/alp/alp2_2021_me1_v_me3_overlaps_pantherGO')
alp2_me1_genes <- read.csv('~/thesis_figs_and_tables/alp/alp2_2021_me1_v_me3_overlaps_FC_UP_overlap.pdf_all_genes_overlaps.csv') %>% dplyr::select(me1)
alp2_me3_genes <- read.csv('~/thesis_figs_and_tables/alp/alp2_2021_me1_v_me3_overlaps_FC_UP_overlap.pdf_all_genes_overlaps.csv') %>% dplyr::select(me3)
panther_go_maker(alp2_me1_genes$me1,'~/thesis_figs_and_tables/alp/alp2_2021_me1_UP_pantherGO')
go_dotplotter('~/thesis_figs_and_tables/alp/alp2_2021_me1_UP_pantherGO')
panther_go_maker(alp2_me3_genes$me3,'~/thesis_figs_and_tables/alp/alp2_2021_me3__UP_pantherGO')
go_dotplotter('~/thesis_figs_and_tables/alp/alp2_2021_me3__UP_pantherGO')



chip_clf_up <- read.table('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff.txt', header=T) %>% filter(FC_KO >1.5, FDR_KO<0.05)

length(chip_clf$FC_KO)
length(chip_clf_up$Chromosome)
trans_clf1$exp <- 'clf28_alp1_rna' 
trans_clf <- filter(trans_clf, gene_id %in% liswt, log2FoldChange>1)
trans_clf <- trans_clf[trans_clf$log2FoldChange >1.5,]

length(trans_clf$gene_id)


liswt <- intersect(chip_clf$feature, trans_clf$gene_id)



ego_rescues <- enrichGO(gene =trans_clf$gene_id, keyType = "TAIR", OrgDb = org.At.tair.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
p <- dotplot(ego_rescues)


trans_alp1 <-read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/alp1-1.csv') %>% dplyr::select(gene_id, log2FoldChange, padj) %>% mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05)
trans_alp1$exp <- 'clf28_alp1_rna'
trans_alp1 <- trans_alp1[trans_alp1$log2FoldChange < -1.5,]
trans_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/alp2-1.csv')%>% dplyr::select(gene_id, log2FoldChange, padj) %>% mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05)
trans_alp2$exp <- 'clf28_alp2_rna'
trans_alp2 <- trans_alp2[trans_alp2$log2FoldChange < -1.5,]
chip_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_transposons.csv') %>% dplyr::select(gene_id, FC_KO, P_KO)%>% dplyr::rename(log2FoldChange=FC_KO,padj=P_KO) %>% dplyr::distinct()%>%filter(padj <0.05)
chip_alp1$exp <- 'clf28_alp1_chip'
chip_alp1 <- chip_alp1[chip_alp1$log2FoldChange > 1.5,]
chip_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_transposons.csv') %>% dplyr::select(gene_id, FC_KO, P_KO)%>% dplyr::rename(log2FoldChange=FC_KO,padj=P_KO)%>% dplyr::distinct() %>% filter(padj <0.05)
chip_alp2$exp <- 'clf28_alp2_chip'
chip_alp2 <- chip_alp2[chip_alp2$log2FoldChange > 1.5,]



5668/length(chip_alp2$gene_id)

trans_alp1 <-read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/alp1-1.csv') %>% dplyr::select(gene_id, log2FoldChange, padj) %>% mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05, gene_id %in% rescue_genelist_clf_down)
trans_alp1$exp <- 'clf28_alp1_rna'
trans_alp1 <- trans_alp1[trans_alp1$log2FoldChange < -1.5,]
trans_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/alp2-1.csv')%>% dplyr::select(gene_id, log2FoldChange, padj) %>% mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05,gene_id %in% rescue_genelist_clf_down)
trans_alp2$exp <- 'clf28_alp2_rna'
trans_alp2 <- trans_alp2[trans_alp2$log2FoldChange < -1.5,]
chip_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_transposons.csv') %>% dplyr::select(gene_id, FC_KO, P_KO)%>% dplyr::rename(log2FoldChange=FC_KO,padj=P_KO) %>% dplyr::distinct()%>%filter(padj <0.05,gene_id %in% rescue_genelist_clf_down)
chip_alp1$exp <- 'clf28_alp1_chip'
chip_alp1 <- chip_alp1[chip_alp1$log2FoldChange > 1.5,]
chip_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_transposons.csv') %>% dplyr::select(gene_id, FC_KO, P_KO)%>% dplyr::rename(log2FoldChange=FC_KO,padj=P_KO)%>% dplyr::distinct() %>% filter(padj <0.05,gene_id %in% rescue_genelist_clf_down)
chip_alp2$exp <- 'clf28_alp2_chip'
chip_alp2 <- chip_alp2[chip_alp2$log2FoldChange > 1.5,]

o <- inner_join(chip_alp2, trans_alp2, by = 'gene_id')
df <- rbind(chip_alp2, trans_alp2) %>% filter(gene_id %in% o$gene_id)
df_for_GO <- df %>% distinct(gene_id,.keep_all = T)

p <- inner_join(chip_alp1, trans_alp1, by = 'gene_id')
df_2 <- rbind(chip_alp1, trans_alp1) %>% filter(gene_id %in% p$gene_id)
df_for_GO <- df %>% distinct(gene_id,.keep_all = T)

q <- rbind(df,df_2)


ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/clf_alp1_alp2_chip_v_rna.pdf',p)


p <- ggplot(q, aes(x=gene_id, y=log2FoldChange, fill=exp)) + geom_bar(stat="identity", color="black", position=position_dodge(), size=0.1) +
    theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"))+
    guides(fill=guide_legend(title="Sample type"))+
    scale_fill_brewer(palette = "Paired")+
    scale_y_continuous(breaks = seq(min(round(df$log2FoldChange)), max(round(df$log2FoldChange)), by = 2))+
    ylab('Fold Change')+
    xlab('Gene name')+
    theme(panel.background = element_blank())


a <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/alp1_alp2_and_rescue_genelist_all_genes_overlaps.csv')


l <- list("AT2G29660","AT4G01920","AT3G44660", "AT4G28490","AT4G37540", "AT5G10945", "AT5G10950", "AT5G15430", "AT2G22980", "AT3G63010")

intersect(l, rescue_genelist_clf_down)

double_bigwigs <- list('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/col_cpm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_cpm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_cpm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2_cpm.bw')
all_peaks <- list('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/col-0_sig_peaks_rescues.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_sig_peaks_rescues.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_alp1_sig_peaks_rescues.csv',
'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_alp2_sig_peaks_rescues.csv')

for(i in 1:length(double_bigwigs)){
    print(i)
    out <- paste0(tools::file_path_sans_ext(all_peaks[[i]]), '_stats')
    print(out)
    get_peak_info(all_peaks[[i]],double_bigwigs[[i]], outname = out)
}
out <- paste0(tools::file_path_sans_ext(all_peaks[[3]]), '_stats')
print(out)
get_peak_info(all_peaks[[3]],double_bigwigs[[3]], outname = out)

clf_down <- read.csv('/Volumes//sesame//ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_both_option.csv') %>% filter(FC_KO < 0.75) %>% filter(feature %in%rescue_genelist_clf_down )
rescue_genelist_clf_down <- read.csv('/Volumes//sesame//ALP_Omics/ChIP//validations/clf_double_mutants/rna_chip_comparison/true_rescues/clf_alp_double_rescue_list.csv')
rescue_genelist_clf_down <- rescue_genelist_clf_down$x


ref6_alp2_genes <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/ref6-alp2_alp2_ref6_as_control_overlaps_all_genes_overlaps.csv') %>% dplyr::select(ref6_alp2)
indices_to_remove <- grep('TE', ref6_alp2_genes$ref6_alp2)
indices_to_remove <- grep(paste0('TE', "|", "-_"), ref6_alp2_genes$ref6_alp2)
ref6_alp2_genes <- ref6_alp2_genes$ref6_alp2[-indices_to_remove]
ref6_alp2_genes <- unique(ref6_alp2_genes)
ref6_alp2_peaks <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_both_fixed.csv') %>% filter(feature %in% true_rescues)



ego_rescues <- enrichGO(gene =l, keyType = "TAIR", OrgDb = org.At.tair.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
p <- dotplot(ego_rescues, showCategory=25)
ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/KEGG_rescue_genes_overlap.pdf',p)
#Shu et al data
swn <- readxl::read_xlsx('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/PLD3-3-e00100-s004_Bound_Genes.xlsx') %>% dplyr::select("SWN unique bound genes") %>% dplyr::rename(genes="SWN unique bound genes")
both <- readxl::read_xlsx('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/PLD3-3-e00100-s004_Bound_Genes.xlsx') %>% dplyr::select("CLF and SWN co-bound genes") %>% dplyr::rename(genes="CLF and SWN co-bound genes")
clf <- readxl::read_xlsx('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/PLD3-3-e00100-s004_Bound_Genes.xlsx') %>% dplyr::select("CLF unique targets") %>% dplyr::rename(genes="CLF unique targets")
clf_peaks <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% dplyr::select(feature, FC_KO) %>% filter(FC_KO > 2 | FC_KO < 0.5)
length(clf_peaks$feature[clf_peaks$FC_KO  > 1])
swn_k_down <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv',skip=1) %>% dplyr::select("swn_K27.decrease")
swn_k_down <- swn_k_down[grep('AT', swn_k_down$swn_K27.decrease),]
clf_k_down <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv',skip=1) %>% dplyr::select("clf_K27.decrease")
clf_k_down <- clf_k_down[grep('AT', clf_k_down$clf_K27.decrease),]


gff <- read.table('/Volumes/sesame/movers/gene_description_20131231(2).txt')

ego_rescues <- enrichGO(gene = rescue_genelist_clf_down, keyType = "TAIR", OrgDb = org.At.tair.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
p <- dotplot(ego_rescues, showCategory=25)
ggsave(plot=p, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/KEGG_enrichment_clf_alp_double_rescue.png')

glist <- list(trans_alp1=trans_alp1$gene_id,trans_alp2=trans_alp2$gene_id,chip_alp1=chip_alp1$gene_id,chip_alp2=chip_alp2$gene_id)
intersection <- find_all_intersections(glist)
rescue_genelist <- unlist(intersection[c("chip_alp1:chip_alp2", "trans_alp1:chip_alp2","trans_alp2:chip_alp1","trans_alp2:chip_alp2","trans_alp1:trans_alp2","trans_alp2:chip_alp1:chip_alp2","trans_alp1:trans_alp2:chip_alp2","trans_alp1:trans_alp2:chip_alp1:chip_alp2")],use.names = F)


my_comparisons <- list( c("clf28", "col-0"), c("clf28", "clf28_alp1"), c("clf28", "clf28_alp2"),
c("clf28_alp1", "clf28_alp2"),c("col-0", "clf28_alp2"),c("col-0", "clf28_alp1"))



#Stats of true rescues
col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/col-0_sig_peaks_rescues_stats') %>% mutate(exp='col-0')
clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_sig_peaks_rescues_stats') %>% mutate(exp='clf28')
clf_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_alp1_sig_peaks_rescues_stats') %>% mutate(exp='clf28_alp1')
clf_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peaks/clf28_alp2_sig_peaks_rescues_stats') %>% mutate(exp='clf28_alp2')
df <- rbind(col,clf,clf_alp1,clf_alp2)

d <- ggviolin(
  df, x = "exp", y = "num_peaks", color = "exp", palette = "Paired",
  add = "jitter", # Add jittered points for better visualization
  fill = "exp"    # Fill the violins with colors based on the "exp" variable
) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 1) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 1) +
  stat_summary(
    fun = "mean",   # Calculate mean for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of mean points
    color = "black" # Set the color of mean points
  ) +
  stat_summary(
    fun = "median", # Calculate median for each "exp" group
    geom = "point", # Add points to the plot
    size = 3,       # Set the size of median points
    color = "red"   # Set the color of median points
  ) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 14, face = "bold")
  )


ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/stats/plots/num_peaks.png',d)


true_rescues <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/clf_alp_double_rescue_list.csv')
true_rescues <- true_rescues$x
alp1_single <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_both_option.csv')%>% dplyr::select(feature, FC_KO, P_KO)%>% dplyr::rename(log2FoldChange=FC_KO,padj=P_KO) %>% dplyr::distinct()%>%filter(padj <0.05,log2FoldChange > 1.5) %>% mutate(exp='alp1')
alp2_single <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_both_option.csv')%>% dplyr::select(feature, FC_KO, P_KO)%>% dplyr::rename(log2FoldChange=FC_KO,padj=P_KO) %>% dplyr::distinct()%>%filter(padj <0.05, log2FoldChange > 1.5) %>% mutate(exp='alp2')



plotgardener::readBigwig(as.character('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2020_ChIP_bams/alp2_rpkm.bw'))



# Load necessary packages
library(rtracklayer)
library(Gviz)

# Define paths to your bigWig files
bigwig_files <- c("sample1.bw", "sample2.bw", "sample3.bw") # Add paths to all your bigWig files

# Create an empty list to store the bigWig objects
bw_list <- list()

# Import the bigWig files into R and store them in the list
for (file in bigwig_files) {
  bw_list[[file]] <- import(file)
}

# Find the maximum value across all samples
max_value <- max(unlist(lapply(bw_list, function(bw) max(values(bw)))))

# Create a function to autoscale a sample
autoscale_sample <- function(sample_bw, max_value) {
  sample_bw$score <- sample_bw$score/max_value
  return(sample_bw)
}

# Apply autoscaling to each sample
scaled_samples <- lapply(bw_list, autoscale_sample, max_value)

# Create a GRanges object for your ChIP-seq peaks
peaks <- import.bed("peaks.bed")

# Create a Gviz plot with autoscaled samples
plotTracks(
  GenomeAxisTrack(),
  AnnotationTrack(range = peaks, name = "ChIP-seq Peaks"),
  lapply(seq_along(scaled_samples), function(i) {
    DataTrack(
      range = scaled_samples[[i]],
      name = paste("Sample", i),
      type = "l",
      ylim = c(0, 1)
    )
  })
  



group_autoscale <- function(sample_dataframes_list) {
  # Create an empty list to store the scaled dataframes
  scaled_dataframes_list <- list()
  
  # Find the maximum value across all samples for each genomic region
  max_values <- lapply(seq_along(sample_dataframes_list[[1]]), function(region_idx) {
    max_val <- max(sapply(sample_dataframes_list, function(sample_df) sample_df[[region_idx]]))
    return(max_val)
  })
  
  # Scale the scores for each sample based on the maximum values
  for (sample_df in sample_dataframes_list) {
    scaled_sample_df <- sample_df
    for (region_idx in seq_along(sample_df)) {
      scaled_sample_df[[region_idx]] <- sample_df[[region_idx]] / max_values[[region_idx]]
    }
    scaled_dataframes_list <- append(scaled_dataframes_list, list(scaled_sample_df))
  }
  
  return(scaled_dataframes_list)
}













opacity_peak_function <- function(coord='gff',namelist=c(),genelist,outprefix,...){
  if(coord == 'gff'){
  gff <- gff <- read.csv('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_GFF3_genes.gff', sep = '\t', header = F) %>%dplyr::filter(V3=='gene')
  }
  else{
    gff <- read.csv(coord)%>% select(chr,feature, start, end) %>% rename('V4'=start, 'V5'=end,'V1'=chr,'V9'=feature) %>% dplyr::distinct()
  }
  bw_list <- list(...)
  df_list <- list()

  
  for(i in 1:length(bw_list)){
    df <- bw_list[i]
    print(as.character(df))
    tracks <- plotgardener::readBigwig(as.character(df))
    tracks$Genotype <- namelist[i]
    df_list[[i]] <- tracks
  }
  max_peak <- 0
  for(gene in genelist){
    print(gene)
    new_df_list <- list()
    final_df_list <- list()
    slim_gff <- gff[grep(gene,gff$V9),]
    print(head(slim_gff))
    start <- slim_gff$V4 -2000
    end <- slim_gff$V5 + 2000
    chr <- slim_gff$V1
    chr <- as.numeric(gsub('Chr','',chr))
    outname <- paste0(outprefix,gene,'.png')
    if (length(start) != 0){
    for(i in 1:length(df_list)){
      df <- as.data.frame(df_list[i])
      subset_df <- df[df$seqnames == slim_gff$V1,]
      seq_start <- which(abs(subset_df$start - start) == min(abs(subset_df$start - start)))
      print(seq_start)
      seq_end <- which(abs(subset_df$start - end) == min(abs(subset_df$start - end)))
      print(seq_end)
      subset_df <- subset_df[seq_start:seq_end,]
      new_df_list[[i]] <- subset_df
    }
    



    for(i in 1:length(new_df_list)){
      df <- as.data.frame(new_df_list[i])
      if (max(df$score) > max_peak){
      max_peak <- max(df$score)
      print(max_peak)
      }
    }
    print(paste('max peak height is ', max_peak))
    #scaled_samples <- lapply(new_df_list, autoscale_sample, max_peak)
    for(i in 1:length(new_df_list)){
      df <- as.data.frame(new_df_list[i])
      print(df$score)
      df$score <-scale(df$score,center=0, scale=max_peak)
      new_df_list[[i]] <- df
    }

    max_peak <- 0
    plot_df <- do.call(rbind, lapply(new_df_list, "[", , c("start", "score","Genotype")))
    plot_df$Genotype <- as.character(plot_df$Genotype)
    p <- ggplot(plot_df, aes(x=start, y=score, fill=Genotype)) + 
    geom_bar(stat="identity", color="transparent", position=position_dodge(), size=0.1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
    ggsave(outname,p)
  }
  }
}

source('/Volumes/sesame/ALP_Omics/ChIP/validations/alp_visualisation_scripts.r')
checknames <- data.frame(Genotype = c('col-0','clf28','clf28 alp1-1', 'clf28 alp2-1'),style = c("plain","italic","italic","italic"))

opacity_peak_function(coord = '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/gff_of_rescues.csv',genelist=rescues$ensembl_gene_id,namelist = checknames,
outprefix = '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/peak_figures/peak_figure_comparison_gff_coords_col_RPGC',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_H3K27me3_RPGC.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_H3K27me3_RPGC.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2_H3K27me3_RPGC.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/col_H3K27me3_RPGC.bw')
clf <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_rpkm.bw'
clf_alp1 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_rpkm.bw'
clf_alp2 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2_rpkm.bw'
col <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/col_rpkm.bw'

gff_for_peaks <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/gff_of_rescues.csv')
glist <- list(clf, clf_alp1, clf_alp2, col)
titles <- list('clf28', 'clf28 alp1', 'clf28 alp2','Col-0')
for (i in 1:length(gff_for_peaks$feature)){
  start <- gff_for_peaks$start[i]
  end <- gff_for_peaks$end[i]
  gene <- gff_for_peaks$feature[i]
  chr <- as.numeric(gsub('Chr','',gff_for_peaks$chr[i]))
  plot_peaks(glist,chr,start-2000,end + 2000,gene, titles, paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/track_figures/', gene,'.png'))
}



###Checking with diffbind and macs2
clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Differential_peaks/res_deseq_2019_clf28_alp2_dba_clf_control_annotated.csv') %>% filter(feature %in% true_rescues)
col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/Differential_peaks/res_deseq_2019_clf28_alp2_dba_col_control_annotated.csv')%>% filter(feature %in% true_rescues)


clf_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_both_option.csv')
up <- filter(clf_alp2, FC_KO > 1)
down <- filter(clf_alp2, FC_KO < 1)


chip_clf_up <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_both_option.csv') %>% filter(FC_KO >2, FDR_KO<0.05)

alp2_single <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_both_option.csv')%>% dplyr::select(feature, FC_KO, P_KO, FDR_KO,peak) %>%filter(feature %in% chip_clf_up$feature)
mean(alp2_single$FC_KO)
alp2_res <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_both_option.csv')%>% dplyr::select(feature, FC_KO, P_KO, FDR_KO,peak) %>% filter(FC_KO > 2)
alp2_diffbind <- read.table('/Volumes/sesame/ALP_Omics/ChIP/Differential_peaks/2020_ChIP/')
mean(alp2_res$FC_KO)

 %>% filter(feature %in% true_rescues)
intersect(alp2_single$feature, both$genes)
ego_rescues <- enrichGO(gene =clf_k_down, keyType = "TAIR", OrgDb = org.At.tair.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
p <- dotplot(ego_rescues, showCategory=25)

diffbind <- read.table('/Volumes/sesame/ALP_Omics/ChIP/Differential_peaks/2019_ChIP/2019_clf28_alp2_dba_ChIP_annotated.txt', header = T)
rescue_diffbind<- diffbind %>% dplyr::filter(feature %in% true_rescues)
intersect(rescue_diffbind$feature, annotated$feature)


alp2_macs3 <- read.csv('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_macs3/clf28_alp2_peaks_annotated.xls')  %>% filter(!(feature %in% true_rescues))
col_macs3 <- read.csv('/Volumes/sesame/joerecovery/Project_folder/alp_omics//2019_macs3/col-0_peaks_annotated.xls')%>% filter(!(feature %in% true_rescues))
alp2_macs3_tip <- read.table('/Volumes/sesame/joerecovery/Project_folder/alp_omics//2019_macs3/clf28_alp2_summits.bed')
col_macs3_tip <- read.table('/Volumes/sesame/joerecovery/Project_folder/alp_omics//2019_macs3/col-0_summits.bed')

# peaks <- alp2_macs3 %>%
#     data.frame() %>%
#     dplyr::rename("chr" = chr, "start" = start, "end" = end)

#   gr <- GenomicRanges::makeGRangesFromDataFrame(
#     peaks,
#     ignore.strand = T, # you sure we need this TRUE?
#     seqnames.field = "chr", # can we get rid of this?
#     start.field = "start",
#     end.field = "end"
#   )

# mart <- biomaRt::useMart(
#     host = "https://plants.ensembl.org", biomart = "plants_mart",
#     dataset = "athaliana_eg_gene"
#   )

# annotated <- ChIPpeakAnno::annotatePeakInBatch(
#     gr,
#     mart,
#     featureType = "TSS",
#     output = "nearestLocation",
#     PeakLocForDistance = "middle"
#   ) %>% Repitools::annoGR2DF() %>% dplyr::select(start,feature)


# out <- dplyr::full_join(annotated, alp2_macs3, by='start')
# write.csv(out, '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_macs3/clf28_alp2_peaks_annotated.xls', quote=F, row.names=F)



  





which(abs(alp2_macs3$start - 10382572) == min(abs(alp2_macs3$start - 10382572)))
which(abs(col_macs3$start - 10382572) == min(abs(col_macs3$start - 10382572)))
#You can't directly use them if experiments have different sequencing depths. You can use '--SPMR' option for 'callpeak' command to get normalized pileup values in million reads. Then these values can be compared.