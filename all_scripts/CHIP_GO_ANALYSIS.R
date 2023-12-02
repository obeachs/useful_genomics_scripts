library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(VennDiagram)
library(RVenn)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(pals)
library(drawProteins)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(yyplot)
library(org.At.tair.db)
library('coriell')


source('~/Salba_RNA/scripts/GO_dotplotter.R')
#Function that takes a genelist (preferably highly significantly up or downregulated genes)
#and makes three barplots 
#if you get an error saying that connection was refused, the number of genes input
# was too big (probably 5000 but hard to find docs)
#Another error can come up if it says curl error, this is likely an internet connection issue
?coriell::panther_go
panther_go_maker <- function(genelist,out_name){
GO_out <- panther_go(gene_list = genelist,organism = 3702,annot_dataset = 'biological_process')
print(head(GO_out))
GG <- GO_out %>% dplyr::select(term,fold_enrichment, fdr) %>% filter(fdr < 0.05)
GG <- GG[!grepl('GO:', GG$term),] %>% arrange(fdr)%>% mutate(term=unlist(term))
write.csv(GG,out_name, quote=F, row.names=F)
}


go_anl <- function(genelist, title='bonk', outprefix = ''){
  annot_list <- c("biological_process", "molecular_function", "cellular_component")
  plot_list <- list()
  for (i in 1:length(annot_list)){
    print(annot_list[i])
    GO_out <- panther_go(gene_list = genelist,organism = 3702,annot_dataset = annot_list[i])
    GO_out <- data.frame(GO_out[!grepl("GO:",GO_out$term),]) %>% filter(fdr < 0.05)
    GO_out <- data.frame(apply(GO_out,2,as.character))
    if (length(GO_out$term) < 1 ) next
    GO_out$fold_enrichment <- as.numeric(GO_out$fold_enrichment)
    if (length(GO_out$term) >20 ){
    GO_out <- head(GO_out,n=20)
    }
    #GO_out <- arrange(GO_out, desc(fold_enrichment)) %>% mutate(name=factor(name,levels=term))
    GO_out$term <- factor(GO_out$term, levels = GO_out$term[order(GO_out$fold_enrichment, decreasing = TRUE)])
    plot_list[[i]] <- ggplot(data = GO_out, aes(x = term, y = fold_enrichment,fill=term)) + geom_col() + geom_hline(yintercept = 1) + 
    geom_abline(intercept = 1, slope = 0) +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
    scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(GO_out$term)))+
     ylab("Fold Enrichment") + 
     ggtitle(paste('Terms for ', annot_list[i], sep=''))
      }
    print(length(plot_list))
    plot <- (ggpubr::ggarrange(plotlist=plot_list, ncol=1,nrow=length(plot_list)))
    if (title!='bonk'){
    plot <- annotate_figure(plot, top = text_grob(title, 
               color = "black", face = "bold", size = 14))
    }
    ggplot2::ggsave(paste(outprefix,'.png',sep=''),plot, height=20, width=30)
  }

#Just doing the biological processes one for this, mostly just to see overlaps
go_anl_list <- function(genelist){
  #if you get an error saying that connection was refused, the number of genes input
  # was too big (probably 5000 but hard to find docs)
  annot_list <- c("biological_process")
  plot_list <- list()
  print(annot_list[1])
  GO_out <- panther_go(gene_list = genelist,organism = 3702,annot_dataset = annot_list[1])
  GO_out <- data.frame(GO_out[!grepl("GO:",GO_out$term),]) %>% filter(fdr < 0.05)
  GO_out <- data.frame(apply(GO_out,2,as.character))
  GO_out$fold_enrichment <- as.numeric(GO_out$fold_enrichment)
  if (length(GO_out$term) >20 ){
  GO_out <- head(GO_out,n=20)
  }
  #GO_out <- arrange(GO_out, desc(fold_enrichment)) %>% mutate(name=factor(name,levels=term))
  
  return(GO_out$term)
}
file <- read.csv('/Users/josephbeegan/Downloads/annotatedchipwithgaplesscalling(1)/kg_transfer/annnotated/AG_3d_IP_annotated_both_option.csv') %>% dplyr::select(feature)
go_anl_list(file)
write.csv(file,'~/Documents/test_go_kg_AP3_3d.txt', quote=F, row.names=F)






source('/Volumes/sesame/joerecovery/scripts/GO_complete.R')
####Making GO heatplots
genelist <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp2-1.csv') %>% dplyr::select(gene_id, padj, log2FoldChange) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5)
clf28_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_midpeak_annotated.csv') %>% dplyr::select(feature, FDR_KO, FC_KO) %>%filter(FDR_KO <0.05) %>% filter(FC_KO >1.5)
clf28_alp2 <- read.csv('/Volumes/sesame/ALP_Omics//ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_midpeak.csv')%>% dplyr::select(feature, FDR_KO, FC_KO)%>%filter(FDR_KO <0.01)%>% filter(FC_KO >1.5)
write.csv(genelist,'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/top_GO/clf_as_control/genelists/clf28_alp1_RNA_genelist.csv', quote=F, row.names=F)

l <- list(clf=clf28$feature,clf_alp1=clf28_alp1$feature,clf_alp2=clf28_alp2$feature)
 
venn.diagram(x = l, filename = '~/Desktop/clf28_doubles_venn.pdf')
genelist <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/top_GO/clf_as_control/genelists/clf28_alp2_RNA_genelist.csv') %>% dplyr::select(gene_id, padj)%>% dplyr::rename('feature'=gene_id, 'FDR_KO'=padj)
get_significant_GO(genelist,gene_col = 'feature',weight_col = 'FDR_KO', save = paste('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/top_GO/clf_as_control/clf28_alp2_RNA__GO',sep = '_'))

files <- list.files('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/top_GO/clf_as_control/genelists', full.names = T, pattern = '.csv')
for (file in files){
  location <- unlist(gregexpr('_genelist', basename(file)))[1]
  outname <-substr(basename(file),1,location)
  outname <- paste(dirname(file),'/',outname,sep = '')
  print(outname)
  if (grepl("RNA", outname)) {
  print("Found 'RNA' in the string.")
  genelist <- read.csv(file) %>% dplyr::select(gene_id,padj) %>% dplyr::filter(padj < 0.05) %>% dplyr::rename('feature'=gene_id, 'FDR_KO'=padj)

  }
  else{
  genelist <- read.csv(file) %>% dplyr::select(feature,FDR_KO) %>% dplyr::filter(FDR_KO < 0.05)
  }
  get_significant_GO(genelist,gene_col = 'feature',weight_col = 'FDR_KO', save = paste(outname,'GO',sep = '_'))
}




gos <- data.frame(GO.ID = c(), Term = c())
for (gene in c("clf28_alp1","clf28_alp2","alp1","alp2","clf28_alp1_RNA", "clf28_alp2_RNA")) {
    dat_tmp <- read.csv(paste0("/Volumes/sesame/ALP_Omics/
    
    hIP/validations/clf_double_mutants/top_GO/clf_as_control/genelists/", gene, '__GOGO_terms.csv'), header = TRUE)
    gos <- rbind.data.frame(gos, dat_tmp[, 1:2])
  }
gos <- unique.data.frame(gos)

gos_df <- data.frame()
for (gene in c("clf28_alp1","clf28_alp2","alp1","alp2","clf28_alp1_RNA", "clf28_alp2_RNA")) {
    dat_tmp <- read.csv(paste0("/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/top_GO/clf_as_control/", gene, '__GOGO_terms.csv'), header = TRUE)
    dat_tmp <- left_join(gos, dat_tmp) %>%
      tibble::add_column(Gene = gene, .before = 1)
    gos_df <- rbind.data.frame(gos_df, dat_tmp)
  }



# Plot significant GOs as heatmaps ---------------------------------------------

# Fisher
gos_df_fisher <- gos_df %>%
  group_by(GO.ID) %>%
  filter(any(Fisher <= 0.01, na.rm = TRUE)) %>%
  ungroup()
gos_df_fisher$Fisher[is.na(gos_df_fisher$Fisher)] <- 1.0
gos_df_fisher$Fisher[gos_df_fisher$Fisher == "<1e-30"] <- 0.0000001
gos_df_fisher$Fisher <- as.numeric(gos_df_fisher$Fisher)

gos_clust <- gos_df_fisher %>%
  dplyr::select(Gene, GO.ID, Fisher) %>%
  tidyr::pivot_wider(names_from = c(Gene), values_from = Fisher) %>%
  tibble::column_to_rownames("GO.ID")
tmp <- data.frame(kmeans(as.matrix(gos_clust), 5, 1000)$cluster) %>%
  tibble::rownames_to_column("GO.ID")
colnames(tmp) <- c("GO.ID", "cluster")


gos_df_fisher2 <- gos_df_fisher %>%
  left_join(tmp) %>%
  arrange(cluster, GO.ID) %>%
  mutate(cluster = factor(cluster)) 
gos_df_fisher2$Fisher[is.na(gos_df_fisher2$Fisher)] <- 1.0
gos_df_fisher2$Fisher <- as.numeric(gos_df_fisher2$Fisher)
gos_df_fisher2$Fisher[is.na(gos_df_fisher2$Fisher)] <- 1.0
gos_df_fisher2$rfisher <- log(gos_df_fisher2$Fisher,10000)




##Makes a heatplot of the GO terms in each of the different genotypes examined
##If there were other timepoints for different genes etc, facet_wrap command 
##in ggplot needs to be used. Change the values in the Fisher column to make the
##opacity of the tiles in geom_tile represent the value correctly - difficult to do 
##if there is a large spread
p <- ggplot(
  mutate(gos_df_fisher2, Term = paste0(stringr::str_sub(Term, end = 50))),
  aes(Gene, reorder(Term, as.numeric(cluster)))
) +
  geom_tile(aes(fill = cluster, alpha = rfisher ), color = "black", size = 0.7) +
  theme_minimal() +
  # scale_fill_gradient(high = "white", low = tropical[1]) +
  scale_fill_manual(values = palette()) +
  scale_alpha(range = c(1, 0)) +
  theme(
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_blank()
  )




# #Use titled lists (i.e genelist <- list(sample_name=sample$feature))
# alp2_ref_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/alp2_me3_diff_ref6_control_annotated_both_fixed.csv') %>% filter(FC_KO > 1.25)
# ref6_alp2_ref_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_both_fixed.csv') %>% filter(FC_KO >1.25)

# alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/alp2_me3_diff_annotated_midpeak_transposons_fixed.csv') %>% filter(FC_KO > 1.25 | FC_KO <0.8)
# ref6_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_alp2_me3_diff_annotated_midpeak_transposons_fixed.csv') %>% filter(FC_KO > 1.25 | FC_KO <0.8)
# ref6 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_me3_diff_annotated_midpeak_transposons_fixed.csv') %>% filter(FC_KO > 1.25 | FC_KO <0.8)


# go_anl(ref6_alp2_ref_control$feature,title = ' ref6 alp2 H3K27me3 DOWN - ref6 control', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/ref6-alp2_ref6_control_DOWN_GO_terms')



# NCBI <- read.table('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_NCBI_GENEID_mapping', sep = '\t', col.names = c('NCBI','TAIR'))
# add_NCBI_ids <- function(results_table){
#   NCBI <- read.table('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_NCBI_GENEID_mapping', sep = '\t', col.names = c('NCBI','feature'))
#   tab <- left_join(results_table, NCBI, by = 'feature')
#   results_table <- tab
#   results_table
# }



# genelist <- list(alp2=alp2$feature[!grepl("TE", alp2$feature)],ref6=ref6$feature[!grepl("TE", ref6$feature)],ref6_alp2=ref6_alp2$feature[!grepl("TE", ref6_alp2$feature)])
# genelist_ref_control <- list(alp2=alp2_ref_control$feature[!grepl("TE", alp2_ref_control$feature)],ref6_alp2=ref6_alp2_ref_control$feature[!grepl("TE", ref6_alp2_ref_control$feature)])
# # Run KEGG analysis
# compKEGG <- compareCluster(gene = genelist_ref_control, 
#                          fun = "enrichKEGG",
#                          organism = "ath",
#                          #OrgDb = "org.At.tair.db",
#                          pvalueCutoff  = 1, 
#                          pAdjustMethod = "BH")


GO_results <- enrichGO(gene = df_alp2$gene_id, 
                       OrgDb = "org.At.tair.db",
                       keyType = 'TAIR',
                       pvalueCutoff = 0.05)


# ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/KEGG_enrichment_comparison.png',dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis alp2 v ref6 alp2 v ref6 (col-0 as control)"))
# ggsave('~/Desktop/testGO.png',dotplot(GO_results, showCategory = 20, title = "KEGG Pathway Enrichment Analysis alp2 v ref6 alp2 v ref6 (col-0 as control)"))

# ##Gene Ontology Analysis

# #Conducting GO analysis and Plotting solely in R.

# #Date: 2023-06-19
# #Load Libs
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# library(tidyverse)
# library(enrichplot)




# #Get the down genes 
# downregulated <- subset(sig, log2FoldChange < 0)

# #Collect the ENSEMBL gene IDs
# genes_to_test <- rownames(downregulated)

# #Get the up genes 
# upregulated <- subset(sig, log2FoldChange > 0)

# #Collect the ENSEMBL gene IDs
# genes_to_test_up <- rownames(upregulated)


# ################################################################################
# #Do the GO analysis (BP)

# str(genes_to_test)
# #Down Genes
# GO_results <- enrichGO(gene = genes_to_test, 
#                        OrgDb = "org.Hs.eg.db",
#                        keyType = "ENSEMBL", ont = "BP",
#                        pvalueCutoff = 0.05)













genelist <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp2-1.csv') %>% dplyr::select(gene_id, padj, log2FoldChange) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5)


genelist <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/clf-28_alp2-1.csv') %>% dplyr::select(gene_id, padj) %>% dplyr::reanm





###Ref6 alp2 
find_all_intersections <- function(list_of_files){
  intersection_total <- gplots::venn(list_of_files,show.plot = FALSE)
  intersection_list <- attributes(intersection_total)$intersections
  intersection_list
}

force_df <- function(listy,maximum){
  c(listy,rep('-',maximum - length(listy)))
}

overlaps_and_venn_output <- function(genes,out_file_prefix,plot_title){
  max_len <- max(unlist(lapply(genes, length)))
  df <- lapply(genes, force_df, max_len)
  write.csv(df,paste(out_file_prefix,'_all_genes.csv', sep = ''), row.names = FALSE, quote = FALSE)
  vennd <- find_all_intersections(genes)
  names(vennd) <- gsub(':', '_' ,names(vennd))
  names(vennd) <-  gsub(' ', '-',names(vennd))
  vennd_max <- max(unlist(lapply(vennd,length)))
  vennd_df <- lapply(vennd, force_df,vennd_max)
  print(names(vennd_df))
  write.csv(vennd_df,paste(out_file_prefix,'_all_genes_overlaps.csv', sep = ''), row.names = FALSE, quote = FALSE)
  p <- ggVennDiagram::ggVennDiagram(genes) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title =plot_title) +  scale_color_brewer(palette = "Paired")+
    theme(
    text = element_text(size = 50),          # Adjust the font size for text
    plot.title = element_text(size = 16),    # Adjust the title font size
    legend.title = element_text(size = 12)  # Adjust the legend title font size
  )+
   guides(fill = 'none')
  ggsave(plot=p,paste(out_file_prefix,'_all_genes_overlaps.pdf', sep = ''), height = 10, width = 10)
}




ref6 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_me3_diff_annotated_midpeak_transposons_fixed.csv') %>% 
filter(FDR_KO < 0.05)%>% filter(FC_KO >2)
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/alp2_me3_annotated_midpeak_transposons_fixed.csv')%>%
filter(FDR < 0.05)%>% filter(log2FoldChange >2)
ref6_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_alp2_me3_annotated_midpeak_transposons_fixed.csv')%>%
filter(FDR < 0.05)%>% filter(log2FoldChange >2)
list <- list(ref6=ref6$feature, alp2=alp2$feature, 'ref6 alp2'=ref6_alp2$feature)
overlaps_and_venn_output(list, '/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps/ref6_alp2_raw_sig_overlaps', 'ref6 alp2 overlaps')



find_all_intersections <- function(list_of_files){
  intersection_total <- gplots::venn(list_of_files,show.plot = FALSE)
  intersection_list <- attributes(intersection_total)$intersections
  intersection_list
}

source('~/Salba_RNA/scripts/GO_dotplotter.R')

ref6_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_both_fixed.csv') %>% 
filter(FC_KO > 1.5)
panther_go_maker(ref6_alp2$feature,'/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/ref6-alp2_ref6_control_UP_pantherGO.csv')
go_dotplotter('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/ref6-alp2_ref6_control_UP_pantherGO.csv')
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/alp2_me3_diff_ref6_control_annotated_both_fixed.csv') %>% 
filter(FC_KO > 2 | FC_KO < 0.5)
rlist <- list('ref6 alp2'=ref6_alp2$feature, alp2=alp2$feature)
overlaps_and_venn_output(rlist, '/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps/ref6_alp2_ref6_control_sig_overlaps', 'ref6 alp2 overlaps - ref6 as control')

ref_control <-read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps/ref6_alp2_ref6_control_sig_overlaps_all_genes_overlaps.csv')
raw <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps/ref6_alp2_raw_sig_overlaps_all_genes_overlaps.csv') %>%
dplyr::select(ref6_alp2_ref6.alp2, ref6_ref6.alp2,ref6.alp2)
raw <- unlist(raw)
raw <- raw[raw !='-']
names(raw) <- NULL
raw <- unique(raw)
eep <- panther_go(raw,organism = 3702,annot_dataset = 'biological_process')
panther_go_maker(raw, '/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/panther_GO/raw_sig_peaks_overlaps.csv')
go_dotplotter('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/panther_GO/raw_sig_peaks_overlaps.csv')

ref_control <-read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/overlaps/ref6_alp2_ref6_control_sig_overlaps_all_genes_overlaps.csv') %>% 
dplyr::select(ref6.alp2_alp2, ref6.alp2)
ref_control <- unlist(ref_control)
ref_control <- ref_control[ref_control !='-']
names(raw) <- NULL
ref_control <- unique(ref_control)
panther_go_maker(ref_control, '/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/panther_GO/ref6_control_peaks_overlaps.csv')
go_dotplotter('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/panther_GO/ref6_control_peaks_overlaps.csv')


panther_go_maker(ref_control$alp2, '/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/panther_GO/ref6_control_alp2_exclusive_peaks_overlaps.csv')
go_dotplotter('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/panther_GO/ref6_control_alp2_exclusive_peaks_overlaps.csv')

alp2 <- alp2 %>% filter(FDR_KO !=1) %>% filter(FC_WT !=1)
clf28_alp_true_rescues <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/clf_alp_double_rescue_list.csv')
panther_go_maker(clf28_alp_true_rescues$x,'/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/clf_alp_double_rescue_list_pantherGO.csv')
go_dotplotter('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/true_rescues/clf_alp_double_rescue_list_pantherGO.csv')

BiocManager::install('yyplot')

##Checking macs3 peaks
source('/Volumes/sesame/ALP_Omics/ChIP/scripts/annotate_peaks_transposable_elements.R')
peaks <- list.files('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_macs3', pattern='annotated', full.names = T)
ref6 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_me3_diff_annotated_midpeak_transposons_fixed.csv') %>% 
filter(FDR_KO < 0.05)%>% filter(FC_KO >2)
ref6 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_me3_diff_annotated_midpeak_transposons_fixed.csv') %>% 
filter(FDR_KO < 0.05)%>% filter(FC_KO > 2)


macs3_ref6 <- read.csv('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_macs3/ref6_peaks.broadPeak_annotated.csv_dedup') %>% 
filter(feature %in% ref6$feature) %>% filter(fc >2)
overlaps_and_venn_output(genes = list(macs3_ref6=macs3_ref6$feature, epic2_ref6=ref6$feature),'/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_sig_raw_peaks_v_macs3_ref6','Significant ref6 raw peaks vs mac3 Ref6 peaks')


ref6_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_alp2_me3_annotated_midpeak_transposons_fixed.csv') %>% 
filter(FDR < 0.05)%>% filter(log2FoldChange >2)
macs3_ref6_alp2 <- read.csv('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_macs3/ref6_alp2_peaks.broadPeak_annotated.csv_dedup') %>% 
filter(feature %in% ref6_alp2$feature)
overlaps_and_venn_output(genes = list(macs3_ref6=macs3_ref6_alp2$feature, epic2_ref6=ref6_alp2$feature),
'/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_sig_diff_peaks_v_macs3_ref6_all_genes',
'Significant ref6 peaks relative to Col-0 vs mac3 ref6 alp2 peaks')


'''Overall, it seems that the epci2 peaks called are all the same as the macs3, just with far more becuase
whatever calculation is used by epic2 must be much more lenient'''


####Finding true alp effectors
ref6_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_both_fixed.csv') %>%
filter(FC_KO < 0.75 | FC_KO > 1.5) %>% distinct()
ref6_alp2_true <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/true_peaks/true_peaks.csv')


alp2_ref_control <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/alp2_me3_diff_ref6_control_annotated_both_fixed.csv') %>%
filter(FC_KO > 1) %>% filter(feature %in% ref6_alp2_true$feature)


'''Only three genes shared between the true ref6 alp2 peaks and the alp2 H3K27me3 increases
relative to ref6'''
ref6_alp2_true_UP <- ref6_alp2_true %>%filter(FC_KO >1)
panther_go_maker(ref6_alp2_true_UP$feature,'/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/true_peaks/H3K27me3_UP_panther_GO')
go_dotplotter('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/true_peaks/H3K27me3_UP_panther_GO')
ref6_alp2_true_DOWN <- ref6_alp2_true %>%filter(FC_KO <1)
panther_go_maker(ref6_alp2_true_DOWN$feature,'/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/true_peaks/H3K27me3_DOWN_panther_GO')
go_dotplotter('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/true_peaks/H3K27me3_DOWN_panther_GO')

rawdata <- data.frame(
  status=c('H3K27me3 UP','H3K27me3 DOWN'),
  value=c(361,216)
)

sigdata <- data.frame(
  status=c('H3K27me3 UP','H3K27me3 DOWN'),
  value=c(102,45)
)

mycols <- c("#FDBF6F","#6FACFC")
raw <- ggplot(rawdata, aes(x = 2, y = value, fill = status)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  labs(title ="ref6 alp2 H3K37me3 relative to ref6",fill = "H3K27me3 Deposition")+
  xlim(0.5, 2.5)+
  annotate(geom = 'text', x = 0.5, y = 0,size=5, label = paste0('Total: ',sum(rawdata$value)))
ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/H3K27me3_donut.pdf', raw, height=10, width=10)

ref6_sig_diff <- ggplot(sigdata, aes(x = 2, y = value, fill = status)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  labs(title ="ref6 alp2 H3K37me3 relative to ref6 - subset of peaks where H3K27me3 increase in ref6 relative to Col-0",fill = "H3K27me3 Deposition")+
  annotate(geom = 'text', x = 0.5, y = 0,size=5, label = paste0('Total: ',sum(sigdata$value)))+
  xlim(0.5, 2.5)
ggsave('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/H3K27me3_ref6_sig_v_col-0_donut.pdf', ref6_sig_diff, height=10, width=10)
