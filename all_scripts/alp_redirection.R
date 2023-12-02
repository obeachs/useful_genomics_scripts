library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gplots)
library(VennDiagram)
library(RVenn)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(pals)
library(clusterProfiler)
library(drawProteins)
library(ggpubr)
library(readxl)
library(knitr)
library(org.At.tair.db)
library('coriell')
library(topGO)
source('~/ChIP_validations.R')

source('/Volumes/sesame/ALP_Omics/ChIP/scripts/annotate_peaks.R')
#Do the alps 'redirect' PRC2 to specific target locations?
#What is the best way to check this?
#Compare the shared genes between the alp single mutants
#and the clf alp double mutants?
#Compare the GO terms of the significant/alp-exclusive genes
#  {
#            "name": "arabidopsis",
#            "taxon_id": 3702,
#            "short_name": "ARATH",
#            "version": "Reference Proteome 2021_03",
#            "long_name": "Arabidopsis thaliana"
#        }

# set a cool color palette
tropical <- c("darkorange", "limegreen", "dodgerblue", "darkorchid", "lightcoral")
palette(tropical)

# GGplot2 theme
theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)

overlaps_and_venn_output <- function(genelist,out_file_prefix,plot_title){
  #Getting the largest amount of features from all of the datasets
  #to be able to make a dataframe from them.
  max_length <- max(unlist(lapply(genelist, length)))
  print(max_length)
  #this requires the force_df function, will fill gaps with NA
  df <- lapply(genelist, force_df, max_length)
  write.csv(df,paste(out_file_prefix,'_all_genes.csv'), row.names = FALSE, quote = FALSE)
  vennd <- find_all_intersections(genelist)
  vennd_max <- max(unlist(lapply(vennd,length)))
  vennd_df <- lapply(vennd, force_df, vennd_max)
  write.csv(vennd_df,paste(out_file_prefix,'_all_genes_overlaps.csv'), row.names = FALSE, quote = FALSE)
  p <- ggVennDiagram::ggVennDiagram(genelist) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title =plot_title)
  ggsave(plot=p,paste(out_file_prefix,'_all_genes_overlaps.pdf'), height = 20, width = 20)
}


get_significant_GO <- function(weighted_genelist,
                               gene_col = "gene_id",
                               weight_col = "pval",
                               ontology = "BP",
                               fun_dplyr::selection = function(x) {
                                 return(x < 0.01)
                               },
                               save = "GO/") {
  
  #formatting the data in a way that is readable for topGO
  #Every column is a gene name and the row of that column is its pvalue
  gene_list <- weighted_genelist[, weight_col]
  names(gene_list) <- weighted_genelist[, gene_col]
  gene_list[is.na(gene_list)] <- 1







  #Creating the topGO object using TAIR
  GOdata <- methods::new(
    "topGOdata",
    description = "BP analysis",
    ontology = ontology,
    allGenes = gene_list,
    nodeSize = 5,
    genedplyr::selectionFun = fun_dplyr::selection,
    annot = annFUN.org,
    mapping = "org.At.tair.db"
  )

#Getting significantly enriched GO terms with Fisher test
  resultFI <- topGO::runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  # resultGL <- topGO::runTest(
  #  GOdata,
  #  algorithm = "weight01",
  #  statistic = "globaltest"
  # )

#Generating a table of the GO term results and writing it to CSV files 
#Taking only significant values
 tabf <- topGO::GenTable(
    GOdata,
    topNodes = max(length(resultFI@score)),
    numChar = 120,
    Fisher = resultFI,
    orderBy = "Fisher",
    ranksOf = "Fisher"
  ) 
  print(head(tabf))
  tabf <- dplyr::filter(tabf,Fisher < 0.05) %>%
    dplyr::filter(Significant > 1) %>%
    dplyr::arrange(Fisher)

    utils::write.csv(tabf, file = paste0(save, "GO_terms.csv"), row.names = F)
    if (length(tabf$Fisher) > 100){
      tabf_short <- dplyr::arrange(tabf, desc(Annotated)) %>% slice(1:100)
      utils::write.csv(tabf_short, file = paste0(save, "GO_terms_top100.csv"), row.names = F)
    }

  # wantedNodes can contain the name of the nodes we want to highlight
  save(GOdata, file = paste0(save, "GOdata.obj"))
   dir.create(paste0(save, "genes_in_GO/"), recursive = TRUE)


  #Making a file showing the GO_terms and the genes associated
  genes <- dplyr::filter(weighted_genelist, FDR_KO < 0.01)$feature
  GO2genes <- topGO::genesInTerm(GOdata)
  GO2genes_subset <- lapply(GO2genes, function(x) x[x %in% genes])
  GO2genes_subset <- GO2genes_subset[lapply(GO2genes_subset, length) > 0]
  tabf <- read.csv(paste0(save, "GO_terms.csv")) %>%
    dplyr::mutate(Term = stringr::str_replace_all(Term, "\\s+", "_")) %>%
    dplyr::mutate(Term = stringr::str_replace_all(Term, "\\/", "-")) %>%
    dplyr::mutate(Term = substr(Term, 1, 40))
  print(head(genes))
  for (go in tabf$GO.ID) {
    go_descr <- dplyr::filter(tabf, GO.ID == go)$Term
    print(go_descr)
    genes <- data.frame(gene_id = GO2genes_subset[[go]]) %>%
      utils::write.csv(file = paste0(save, "genes_in_GO/", go, "_", go_descr, ".csv"), row.names = FALSE)
    print(genes)
    print(paste0(save, "genes_in_GO/", go, "_", go_descr, ".csv"))
  #This generates the ugly plots/flowcharts
  #topGO::printGraph(GOdata, resultFI, firstSigNodes = 5, useInfo = "all", fn.prefix = paste0(save, "FI_flow_chart_first5"), pdfSW = TRUE)
  }
}




clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_both_option.csv') %>% filter(FC_KO < 0.8)
clf_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp1_me3_diff_annotated_both_option.csv') %>% filter(FC_KO < 0.8)
clf_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_alp2_me3_diff_annotated_both_option.csv') %>% filter(FC_KO < 0.8)
alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_both_option.csv') %>% filter(FC_KO < 0.8)
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_both_option.csv') %>% filter(FC_KO < 0.8)
all <- list("clf"=clf$feature,"clf_alp1"=clf_alp1$feature,"clf_alp2"=clf_alp2$feature,"alp1"=alp1$feature,"alp2"=alp2$feature)



counts <- counts[counts$geneid %in% list,]

force_df <- function(listy,maximum){
  c(listy,rep('-',maximum - length(listy)))
}

find_all_intersections <- function(list_of_files){
  intersection_total <- gplots::venn(list_of_files,show.plot = FALSE)
  intersection_list <- attributes(intersection_total)$intersections
  intersection_list
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
  annot_list <- c("biological_process", "molecular_function", "cellular_component")
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
overlaps_and_venn_output <- function(genes,out_file_prefix,plot_title){
  max_len <- max(unlist(lapply(genes, length)))
  df <- lapply(genes, force_df, max_len)
  write.csv(df,paste(out_file_prefix,'_all_genes.csv', sep = ''), row.names = FALSE, quote = FALSE)
  vennd <- find_all_intersections(genes)
  vennd_max <- max(unlist(lapply(vennd,length)))
  vennd_df <- lapply(vennd, force_df,vennd_max)
  write.csv(vennd_df,paste(out_file_prefix,'_all_genes_overlaps.csv', sep = ''), row.names = FALSE, quote = FALSE)
  p <- ggVennDiagram::ggVennDiagram(genes) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title =plot_title)
  ggsave(plot=p,paste(out_file_prefix,'_all_genes_overlaps.pdf', sep = ''), height = 20, width = 20)
}

GO_out <- panther_go(gene_list = alp1_me3_up$feature,organism = 3702,annot_dataset ="biological_process")
all_match <- find_all_intersections(all)
length(all_match$'clf:clf_alp1:clf_alp2:alp1:alp2')
genes <- c(all_match$'alp1')
go_anl(alp1)
pal::pal(colorRampPalette(brewer.pal(20, "Paired")))



#Rescued genes by alp double mutants
list.files('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/', pattern = d)

alp1_me3_up <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_annotated_both.csv') %>% filter(FC_KO > 1.2)
alp1_me3_up_go <-go_anl_list(alp1_me3_up$feature)
go_anl(alp1_me3_up$feature,title='clf alp1 H3K27me3 up 1.2x FC', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/GO_terms/clf_alp1_H3K27me3_up_GO')
alp1_me3_down <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_annotated_both.csv') %>% filter(FC_KO < 0.8)
alp1_me3_down_go <-go_anl_list(alp1_me3_down$feature)
go_anl(alp1_me3_down$feature,title='clf alp1 H3K27me3 down 1.2x FC', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/GO_terms/clf_alp1_H3K27me3_down_GO')
alp2_me3_up <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_annotated_both.csv') %>% filter(FC_KO > 1.3)
alp2_me3_up_go <-go_anl_list(alp2_me3_up$feature)
go_anl(alp2_me3_up$feature,title='clf alp2 H3K27me3 up 1.2x FC', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/GO_terms/clf_alp2_H3K27me3_up_GO')
alp2_me3_down <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_annotated_both.csv') %>% filter(FC_KO < 0.75)
alp2_me3_down_go <-go_anl_list(alp2_me3_down$feature)
go_anl(alp2_me3_down$feature,title='clf alp1 H3K27me3 down 1.2x FC', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/GO_terms/clf_alp2_H3K27me3_down_GO')

GO <- list(alp1_up=alp1_me3_up_go,alp2_up=alp2_me3_up_go, alp1_down=alp1_me3_down_go, alp2_down=alp2_me3_down_go)
GO<- find_all_intersections(GO)


panther_go(gene_list = alp2_me3_up$feature,organism = 3702,annot_dataset = 'biological_process')


salp1_up <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_both_option.csv')%>% filter(FC_KO >1.2)
salp1_down <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_both_option.csv')%>% filter(FC_KO < 0.8)
salp2_up <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_both_option.csv')%>% filter(FC_KO > 1.2)
salp2_down <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_both_option.csv')%>% filter(FC_KO < 0.8)
go_anl(salp1_up$feature, title = 'alp1 H3K27me3 up 1.2x FC', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/GO_terms/alp1_me3_up_GO')
go_anl(salp2_up$feature, title = 'alp2 H3K27me3 up 1.2x FC', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/GO_terms/alp2_me3_up_GO')
go_anl(salp1_down$feature, title = 'alp1 H3K27me3 down 0.8x FC', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/GO_terms/alp1_me3_down_GO')
go_anl(salp2_down$feature, title = 'alp2 H3K27me3 down 0.8x FC', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/GO_terms/alp2_me3_down_GO')

sGO <- list(salp1_up= salp1_up_go,salp2_up= salp2_up_go,salp1_down=salp1_down_go,salp2_down=salp2_down_go)
sGO <- find_all_intersections(sGO)

out_file_prefix,
                           chr_col,
                           start_col,
                           end_col,

annotate_peaks(peaks_file='/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/col-0_me3.txt',
out_file_prefix='/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/col-0_me3_annotated_both',
start_col='Start', end_col='End', chr_col='Chromosome')







#TEs
clf <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/transposon_check/tair_transposable_elements/clf28_me3_diff_annotated_tair_TE_midpeak.csv')
clf_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/transposon_check/tair_transposable_elements/clf28_alp1_me3_diff_annotated_tair_TE_midpeak.csv')
clf_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/transposon_check/tair_transposable_elements/clf28_alp2_me3_diff_annotated_tair_TE_midpeak.csv')
alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/transposon_check/tair_transposable_elements/alp1_me3_diff_annotated_tair_TE_midpeak.csv')
alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/transposon_check/tair_transposable_elements/alp2_me3_diff_annotated_tair_TE_midpeak.csv')
col <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/col-0_me3_annotated_both.csv')

clf_up <- filter(clf, FC_KO > 1.5)
100*(length(clf_up$feature)/length(clf_up$feature[grep('TE', clf_up$feature)]))
clf_alp1_up <- filter(clf_alp1, FC_KO > 1.5)
100*(length(clf_alp1_up$feature[grep('TE', clf_alp1_up$feature)])/length(clf_alp1_up$feature))
clf2_alp1_up <- filter(clf_alp2, FC_KO > 1.5)
100*(length(clf_up$feature)/length(clf_up$feature[grep('TE', clf_up$feature)]))
alp1_up <- filter(alp1, FC_KO > 1.5)
alp2_up <- filter(alp2, FC_KO > 1.5)
clf_down <- filter(clf, FC_KO <0.75)
clf_alp1_down <- filter(clf_alp1, FC_KO <0.75)
100*(length(clf_alp1_down$feature[grep('TE', clf_alp1_down$feature)])/length(clf_alp1_down$feature))
clf2_alp1_down <- filter(clf_alp2, FC_KO <0.75)
alp1_down <- filter(alp1, FC_KO <0.75)
alp2_down <- filter(alp2, FC_KO <0.75)

100*(length(col$feature[grep('TE', col$feature)])/length(col$feature))
col

 
#Looking at the overlap between the clf28 alp double mutants
#ChIP and RNAseq - are the H3K27me3 changes changing the expression
#of the genes?

"""For this section only looking at the clf28 as control"""
clf28_alp1_rna <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/alp1-1.csv') %>% dplyr::dplyr::select(gene_id, log2FoldChange, padj) %>% mutate(gene_id = as.character(gene_id))
clf28_alp1_rna$exp <- 'clf28_alp1_rna'
clf28_alp1_rna <- clf28_alp1_rna[(clf28_alp1_rna$log2FoldChange > 1.5 | clf28_alp1_rna$log2FoldChange < -1.5),]
clf28_alp2_rna <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/alp2-1.csv')%>% dplyr::dplyr::select(gene_id, log2FoldChange, padj) %>% mutate(gene_id = as.character(gene_id))
clf28_alp2_rna$exp <- 'clf28_alp2_rna'
clf28_alp2_rna <- clf28_alp2_rna[(clf28_alp2_rna$log2FoldChange > 1.5 | clf28_alp2_rna$log2FoldChange < -1.5),]


clf28_alp1_chip <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_transposons.csv') %>% dplyr::dplyr::select(gene_id, FC_KO, P_KO)%>% dplyr::dplyr::rename(log2FoldChange=FC_KO,padj=P_KO)
clf28_alp1_chip$exp <- 'clf28_alp1_chip'
clf28_alp1_chip <- clf28_alp1_chip[(clf28_alp1_chip$log2FoldChange > 1.5 | clf28_alp1_chip$log2FoldChange < 0.75),] %>% mutate(log2FoldChange=log(log2FoldChange,2))
clf28_alp2_chip <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_transposons.csv') %>% dplyr::dplyr::select(gene_id, FC_KO, P_KO)%>% dplyr::dplyr::rename(log2FoldChange=FC_KO,padj=P_KO)
clf28_alp2_chip$exp <- 'clf28_alp2_chip'
clf28_alp2_chip <- clf28_alp2_chip[(clf28_alp2_chip$log2FoldChange > 1.5 | clf28_alp2_chip$log2FoldChange < 0.5),] %>% mutate(log2FoldChange=log(log2FoldChange,2))

alp1_df <- inner_join(clf28_alp1_rna, clf28_alp1_chip, by = 'gene_id')
df <- rbind(clf28_alp1_rna,clf28_alp1_chip)%>% filter(gene_id %in% alp1_df$gene_id)
alp2_df <- inner_join(clf28_alp2_rna, clf28_alp2_chip, by = 'gene_id')
df_alp2 <- rbind(clf28_alp2_rna,clf28_alp2_chip) %>% filter(gene_id %in% alp2_df$gene_id)
df_1 <- rbind(clf28_alp1_rna,clf28_alp1_chip)
df_2 <- rbind(clf28_alp2_rna,clf28_alp2_chip)
df_2 <- inner_join(clf28_alp2_rna,clf28_alp2_chip, by = 'gene_id')
full <- rbind(df,df_alp2)
full <- rbind(df_1,df_2)

full$padj <- as.numeric(full$padj)
full <- full[full$padj <1,]

glist <- list(trans_alp1=trans_alp1$gene_id,trans_alp2=trans_alp2$gene_id,chip_alp1=chip_alp1$gene_id,chip_alp2=chip_alp2$gene_id)
intersection <- find_all_intersections(glist)

#ONLY TRANSCRIPTIONAL AND H3K27me3 RESCUES]
trans_alp1 <-read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/alp1-1.csv') %>% dplyr::dplyr::select(gene_id, log2FoldChange, padj) %>% mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05)
trans_alp1$exp <- 'clf28_alp1_rna'
trans_alp1 <- trans_alp1[trans_alp1$log2FoldChange < -1.5,]
trans_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/RNASeq/DE_results/Bennett/noshrink/reference_clf-28/alp2-1.csv')%>% dplyr::dplyr::select(gene_id, log2FoldChange, padj) %>% mutate(gene_id = as.character(gene_id))%>% filter(padj <0.05)
trans_alp2$exp <- 'clf28_alp2_rna'
trans_alp2 <- trans_alp2[trans_alp2$log2FoldChange < -1.5,]
chip_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_transposons.csv') %>% dplyr::dplyr::select(gene_id, FC_KO, P_KO)%>% dplyr::dplyr::rename(log2FoldChange=FC_KO,padj=P_KO) %>% dplyr::distinct()%>%filter(padj <0.05)
chip_alp1$exp <- 'clf28_alp1_chip'
chip_alp1 <- chip_alp1[chip_alp1$log2FoldChange > 1.5,]
chip_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_transposons.csv') %>% dplyr::dplyr::select(gene_id, FC_KO, P_KO)%>% dplyr::dplyr::rename(log2FoldChange=FC_KO,padj=P_KO)%>% dplyr::distinct() %>% filter(padj <0.05)
chip_alp2$exp <- 'clf28_alp2_chip'
chip_alp2 <- chip_alp2[chip_alp2$log2FoldChange > 1.5,]
o <- inner_join(chip_alp2, trans_alp2, by = 'gene_id')
df <- rbind(chip_alp2, trans_alp2) %>% filter(gene_id %in% o$gene_id)
df_for_GO <- df %>% distinct(gene_id,.keep_all = T)




glist <- list(trans_alp1=trans_alp1$gene_id,trans_alp2=trans_alp2$gene_id,chip_alp1=chip_alp1$gene_id,chip_alp2=chip_alp2$gene_id)
intersection <- find_all_intersections(glist)
rescue_genelist <- unlist(intersection[c("chip_alp1:chip_alp2", "trans_alp1:chip_alp2","trans_alp2:chip_alp1","trans_alp2:chip_alp2","trans_alp1:trans_alp2","trans_alp2:chip_alp1:chip_alp2","trans_alp1:trans_alp2:chip_alp2","trans_alp1:trans_alp2:chip_alp1:chip_alp2")],use.names = F)
#Check for overlaps with SWN bound regions
swn <- readxl::read_xlsx('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/PLD3-3-e00100-s004_Bound_Genes.xlsx') %>% dplyr::dplyr::select("SWN unique bound genes") %>% dplyr::dplyr::rename(genes="SWN unique bound genes")
both <- readxl::read_xlsx('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/PLD3-3-e00100-s004_Bound_Genes.xlsx') %>% dplyr::dplyr::select("CLF and SWN co-bound genes") %>% dplyr::dplyr::rename(genes="CLF and SWN co-bound genes")
clf <- readxl::read_xlsx('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/PLD3-3-e00100-s004_Bound_Genes.xlsx') %>% dplyr::dplyr::select("CLF unique targets") %>% dplyr::dplyr::rename(genes="CLF unique targets")
clf_peaks <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/default_params/clf28_me3_diff_annotated_midpeak.csv') %>% dplyr::dplyr::select(feature, FC_KO) %>% filter(FC_KO > 2 | FC_KO < 0.5)
length(clf_peaks$feature[clf_peaks$FC_KO  > 1])
swn_k_down <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv',skip=1) %>% dplyr::dplyr::select("swn_K27.decrease")
swn_k_down <- swn_k_down[grep('AT', swn_k_down$swn_K27.decrease),]
clf_k_down <- read.csv('/Volumes/sesame/ALP_Omics/ext_datasets/Shu2019_CLF_SWN_targets/CLF-SWN_H3K27me3.csv',skip=1) %>% dplyr::dplyr::select("clf_K27.decrease")
clf_k_down <- clf_k_down[grep('AT', clf_k_down$clf_K27.decrease),]
length(intersect(rescue_genelist, clf_k_down))


rescue_genelist[order(names(setNames(rescue_genelist, rescue_genelist)))]




###DOES this track with the alp single mutants?
alp1_single <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp1_me3_diff_annotated_both_option.csv')%>% dplyr::dplyr::select(feature, FC_KO, P_KO)%>% dplyr::dplyr::rename(log2FoldChange=FC_KO,padj=P_KO) %>% dplyr::distinct()%>%filter(padj <0.05,log2FoldChange > 1.5) %>% mutate(exp='alp1')
alp2_single <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2020_ChIP/epic2/default_params/alp2_me3_diff_annotated_both_option.csv')%>% dplyr::dplyr::select(feature, FC_KO, P_KO)%>% dplyr::dplyr::rename(log2FoldChange=FC_KO,padj=P_KO) %>% dplyr::distinct()%>%filter(padj <0.05, log2FoldChange > 1.5) %>% mutate(exp='alp2')
intersect(singles$feature, rescue_genelist)

genelist,out_file_prefix,plot_title
singles <- rbind(alp1_single, alp2_single)
for_venn <- list(rescues=rescue_genelist, alp1=alp1_single$feature, alp2=alp2_single$feature)
overlaps_and_venn_output(for_venn, out_file_prefix = '/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/alp1_alp2_and_rescue_genelist',plot_title = 'Overlap of ALP single mutant H3K27me3 Peaks and the Rescue genelist')

intersect(ego@result$Description,ego_singles@result$Description )
singles_sig <- filter(singles, log2FoldChange >1.5)
singles_go <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/alp1_alp2_rescue_genelist_overlap_topGOGO_terms.csv')
ego_singles <- enrichGO(gene = unique(singles_sig$feature), 
                    keyType = "TAIR", 
                    OrgDb = org.At.tair.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

p <- dotplot(ego_singles, showCategory=5000)
ego_rescues <- enrichGO(gene = rescue_genelist, 
                    keyType = "TAIR", 
                    OrgDb = org.At.tair.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

p <- dotplot(ego, showCategory=5000)
ggsave(plot = p, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/RNA_ChIP_merge_KEGG_BP_enrich_new.png')
p <- ggplot(singles_sig, aes(x=feature, y=log2FoldChange, fill=exp)) + geom_bar(stat="identity", color="black", position=position_dodge(), size=0.1) +
    theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"))+
    guides(fill=guide_legend(title="Sample type"))+
    scale_fill_brewer(palette = "Paired")+
    scale_y_continuous(breaks = seq(min(round(df$log2FoldChange)), max(round(df$log2FoldChange)), by = 2))+
    ylab('Fold Change')+
    xlab('Gene name')+
    theme(panel.background = element_blank())
print(p)
get_significant_GO(dplyr::dplyr::select(singles, feature, padj) %>% dplyr::rename(FDR_KO=padj),gene_col = 'feature', weight_col = 'FDR_KO', save = '/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/alp1_alp2_rescue_genelist_overlap_topGO')


go_anl(genelist = df$gene_id, title = 'clf28 alp double mutants using clf28 as control RNA and ChIP', outprefix = '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/RNA_ChIP_clf-doubles_clf_as_control_GO')
get_significant_GO(dplyr::dplyr::select(df_for_GO, gene_id, padj) %>% dplyr::dplyr::rename(FDR_KO=padj, feature=gene_id),gene_col = 'feature', weight_col = 'FDR_KO', save = '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rescue_genelist_topGO')

ego <- enrichGO(gene = rescue_genelist, 
                    keyType = "TAIR", 
                    OrgDb = org.At.tair.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

p <- dotplot(ego, showCategory=50)
ggsave(plot = p, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/RNA_ChIP_merge_KEGG_BP_enrich_new.png')

p <- ggplot(df, aes(x=gene_id, y=log2FoldChange, fill=exp)) + geom_bar(stat="identity", color="black", position=position_dodge(), size=0.1) +
    theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"))+
    guides(fill=guide_legend(title="Sample type"))+
    scale_fill_brewer(palette = "Paired")+
    scale_y_continuous(breaks = seq(min(round(df$log2FoldChange)), max(round(df$log2FoldChange)), by = 2))+
    ylab('Fold Change')+
    xlab('Gene name')+
    theme(panel.background = element_blank())
print(p)

ggsave(plot = p, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/RNA_ChIP_clf-doubles_clf_as_control_barplot.pdf')

p <- ggplot(df, aes(exp,gene_id)
) +
  #facet_wrap(~exp, ncol = 2) +
  geom_tile(aes(fill = exp, alpha =log2FoldChange), color = "black", size = 0.1) +
  theme_minimal() +
  # scale_fill_gradient(high = "white", low = tropical[1]) +
  scale_fill_manual(values = palette()) +
  scale_alpha(range = c(1, 0)) +
     xlab('Experiment') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(plot = p, '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/clf28_alp1_RNA_ChIP_clf_as_control_heatplot.pdf')



gff <- read.csv('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_GFF3_genes.gff', sep = '\t', header = F) %>% dplyr::filter(V3=='gene')
#Looking at the difference in the peaks that show expression changes compared to clf28
for (name in full$gene_id){
  df <- gff[grep(name,gff$V9),]
  start <- df$V4
  end <- df$V5
  chr <- df$V1
  chr <- as.numeric(gsub('Chr','',chr))
  print(start)
  print(end)
  max_yo <- 0
  outname <- paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/',name,'_significant_tracks_bamcompare.png')
  if (is.integer(start)){
    bigwig_list <- list('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28.bw','/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1-1.bw',
    '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1-1.bw')
    print(max_height)
    for(wig in bigwig_list){
    tracks <- plotgardener::readBigwig(wig)
    tracks$seqnames <- gsub('Chr','', tracks$seqnames)
    chrom_name <- chr
    subset_df <- tracks[tracks$seqnames == chr,]
    seq_start <- which(abs(subset_df$start - start) == min(abs(subset_df$start - start)))
    seq_end <- which(abs(subset_df$start - end) == min(abs(subset_df$start - end)))
    subset_df <- subset_df[seq_start:seq_end,]
    print(max(subset_df$score))
    new_y <- max(subset_df$score)
    if (as.numeric(max(subset_df$score)) > max_yo){
      max_yo <- max(subset_df$score)
    }
    }
    print(max_yo)
    clf28 <- plot_peaks(bigwig = '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28.bw',chrom = chr,start = start-1000,end = end+1000, genename = 'clf', max_y = max_yo)
    clf28_alp1_1 <- plot_peaks(bigwig = '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1-1.bw',chrom = 4,start = start-1000,end = end+1000,genename = 'clf alp1-1', max_y = max_yo)
    clf28_alp2_1 <- plot_peaks(bigwig = '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2-1.bw',chrom = 4,start = start-1000,end = end+1000,genename = 'clf alp2-1', max_y = max_yo)
    p <- ggpubr::ggarrange(clf28,clf28_alp1_1,clf28_alp2_1,nrow = 3,ncol =1 )
    ggsave(outname, p)
  } # Apply is.integer function # FALSE.
}

for (name in unique(full$gene_id)){
  df <- gff[grep(name,gff$V9),]
  start <- df$V4
  end <- df$V5
  chr <- df$V1
  chr <- as.numeric(gsub('Chr','',chr))
  print(start)
  print(chr)
  print(end)
  max_yo <- 0
  tracks <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_rpkm.bw')
  atracks <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_rpkm.bw')
  outname <- paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/',name,'_significant_tracks_bamcompare_opacity_rpkm_alp1.png')
  if (start > 1){
    tracks$seqnames <- gsub('Chr','', tracks$seqnames)
    chrom_name <- chr
    clfsubset_df <- tracks[tracks$seqnames == chr,]
    seq_start <- which(abs(clfsubset_df$start - start) == min(abs(clfsubset_df$start - start)))
    seq_end <- which(abs(clfsubset_df$start - end) == min(abs(clfsubset_df$start - end)))
    clfsubset_df <- clfsubset_df[seq_start:seq_end,]
    atracks$seqnames <- gsub('Chr','', atracks$seqnames)
    chrom_name <- chr
    alp2subset_df <- atracks[atracks$seqnames == chr,]
    seq_start <- which(abs(alp2subset_df$start - start) == min(abs(alp2subset_df$start - start)))
    seq_end <- which(abs(alp2subset_df$start - end) == min(abs(alp2subset_df$start - end)))
    alp2subset_df <- alp2subset_df[seq_start:seq_end,]
    clfsubset_df$exp <- 'clf'
    alp2subset_df$exp <- 'clf-alp1'
    final_df <- rbind(clfsubset_df, alp2subset_df)
    print(head(final_df))
    p <- ggplot(final_df, aes(x=start, y=score, fill=exp)) + 
    geom_bar(stat="identity", color="transparent", position=position_dodge(), size=0.1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
    ggsave(outname, p)
  }
}

for (name in unique(full$gene_id)){
  df <- gff[grep(name,gff$V9),]
  start <- df$V4
  end <- df$V5
  chr <- df$V1
  chr <- as.numeric(gsub('Chr','',chr))
  print(start)
  print(chr)
  print(end)
  max_yo <- 0
  tracks <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_rpkm.bw')
  atracks <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2_rpkm.bw')
  btracks <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_rpkm.bw')
  outname <- paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/',name,'_significant_tracks_bamcompare_opacity_rpkm_alp1_alp2.png')
  if (start != integer(0)){
    chrom_name <- chr
    
    tracks$seqnames <- gsub('Chr','', tracks$seqnames)
    clfsubset_df <- tracks[tracks$seqnames == chr,]
    seq_start <- which(abs(clfsubset_df$start - start) == min(abs(clfsubset_df$start - start)))
    seq_end <- which(abs(clfsubset_df$start - end) == min(abs(clfsubset_df$start - end)))
    clfsubset_df <- clfsubset_df[seq_start:seq_end,]

    atracks$seqnames <- gsub('Chr','', atracks$seqnames)
    alp2subset_df <- atracks[atracks$seqnames == chr,]
    seq_start <- which(abs(alp2subset_df$start - start) == min(abs(alp2subset_df$start - start)))
    seq_end <- which(abs(alp2subset_df$start - end) == min(abs(alp2subset_df$start - end)))
    alp2subset_df <- alp2subset_df[seq_start:seq_end,]
    tracks$seqnames <- gsub('Chr','', tracks$seqnames)

    btracks$seqnames <- gsub('Chr','', btracks$seqnames)
    bclfsubset_df <- btracks[btracks$seqnames == chr,]
    seq_start <- which(abs(bclfsubset_df$start - start) == min(abs(bclfsubset_df$start - start)))
    seq_end <- which(abs(bclfsubset_df$start - end) == min(abs(bclfsubset_df$start - end)))
    bclfsubset_df <- bclfsubset_df[seq_start:seq_end,]

    clfsubset_df$exp <- 'clf'
    alp2subset_df$exp <- 'clf-alp2'
    bclfsubset_df$exp <- 'clf-alp1'
    final_df <- rbind(clfsubset_df, alp2subset_df, bclfsubset_df)
    print(head(final_df))
    p <- ggplot(final_df, aes(x=start, y=score, fill=exp)) + 
    geom_bar(stat="identity", color="transparent", position=position_dodge(), size=0.1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
    ggsave(outname, p)
  }
}


#dplyr::selecting the top 30 genes in terms of FC from the ref6-alp2 using ref6 as control to see what the alps are potentially involved in.
ref6_alp2_genes <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_as_control/ref6_alp2_me3_diff_ref6_control_annotated_both.csv') %>% filter(FC_KO > 1.5)%>% dplyr::arrange(desc(FC_KO)) %>% dplyr::dplyr::select(feature) %>% dplyr::distinct(feature) %>% slice(1:30)
alp1_alp2_genes <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/alp1_alp2/differential_peaks/alp2_diff_sig_peaks.csv')
#alp2_genes <- unique(read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_me3_diff_annotated_midpeak_transposons_fixed.csv') %>%

genes <- intersect(ref6$feature, ref6_alp2_genes$feature)
opacity_peak_function <- function(coord='gff',namelist=c(),genelist,outprefix,...){
  if(coord == 'gff'){
  gff <- gff <- read.csv('/Volumes/sesame/joerecovery/genomes/TAIR/TAIR10_GFF3_genes.gff', sep = '\t', header = F) %>%dplyr::filter(V3=='gene')
  }
  else{
    gff <- read.csv(coord)%>% dplyr::select(chr,feature, start, end) %>% dplyr::rename('V4'=start, 'V5'=end,'V1'=chr,'V9'=feature) %>% dplyr::distinct()
  }
  bw_list <- list(...)
  df_list <- list()
  for(i in 1:length(bw_list)){
    df <- bw_list[i]
    print(as.character(df))
    tracks <- plotgardener::readBigwig(as.character(df))
    tracks$Genotype <- namelist$Genotype[i]
    tracks$style <- namelist$style[i]
    df_list[[i]] <- tracks
  }
  for(gene in genelist){
    new_df_list <- list()
    slim_gff <- gff[grep(gene,gff$V9),]
    print(head(slim_gff))
    start <- slim_gff$V4 -5000
    print(start)
    end <- slim_gff$V5 + 5000
    print(end)
    chr <- slim_gff$V1
    chr <- as.numeric(gsub('Chr','',chr))
    print(chr)
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
    plot_df <- do.call(rbind, lapply(new_df_list, "[", , c("start", "score","Genotype", "style")))
    plot_df$Genotype <- as.character(plot_df$Genotype)
    print(plot_df)
    p <- ggplot(plot_df, aes(x=start, y=score, fill=Genotype)) + 
    geom_bar(stat="identity", color="transparent", position=position_dodge(), size=0.1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    
    scale_fill_brewer(palette ='Set2')
    ggsave(outname,p)
  }
  }
}

checknames <- data.frame(Genotype = c('ref6','ref6 alp2-1','alp2-1', 'Col-0'),
                            style = c("italic","italic","italic","plain"))

r <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_alp2_me3_diff_annotated_midpeak_transposons_fixed_ref_2xUP.csv')

#Ref-6 alp2-1 with ref-6 as control peaks
opacity_peak_function(coord='/Volumes/sesame/ALP_Omics/ChIP/2021_ChIP/epic2/ref6_alp2_me3_diff_annotated_midpeak_transposons_fixed_ref_2xUP.csv',namelist=checknames,outprefix='/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/col-0_control/peak_figures_ref6_up_v_col/peak_figure_comparison_peak_coords_',genelist=alp1_alp2_genes$feature,'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_rpkm.bw',
'/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_alp2_rpkm.bw','/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/alp2_rpkm.bw','/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/col_rpkm.bw')
  #clf alp double mutants
opacity_peak_function(coord = 'gff',genelist= full$gene_id,namelist = list('clf28','clf28 alp1-1','clf28 alp2-1', 'col-0'), outprefix = '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/peak_figure_comparison_gff_coords_col','/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_rpkm.bw','/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1_rpkm.bw','/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2_rpkm.bw','/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/col_rpkm.bw')


tracks <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28.bw')
tracks$seqnames <- gsub('Chr','', tracks$seqnames)
chrom_name <- 1
subset_df <- tracks[tracks$seqnames == 1,]
seq_start <- which(abs(subset_df$start - 29000527) == min(abs(subset_df$start - 29000527)))
seq_end <- which(abs(subset_df$start - 29008920) == min(abs(subset_df$start - 29008920)))
subset_df <- subset_df[seq_start:seq_end,]

atracks <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp2-1.bw')
atracks$seqnames <- gsub('Chr','', atracks$seqnames)
chrom_name <- 1
asubset_df <- atracks[atracks$seqnames == 1,]
seq_start <- which(abs(asubset_df$start - 29000527) == min(abs(asubset_df$start - 29000527)))
seq_end <- which(abs(asubset_df$start - 29008920) == min(abs(asubset_df$start - 29008920)))
asubset_df <- asubset_df[seq_start:seq_end,]
subset_df$exp <- 'clf'
asubset_df$exp <- 'clf-alp2'
doub <- rbind(subset_df, asubset_df)

btracks <- plotgardener::readBigwig('/Volumes/sesame/joerecovery/Project_folder/alp_omics/2019_ChIP_bams/clf28_alp1-1.bw')
btracks$seqnames <- gsub('Chr','', btracks$seqnames)
chrom_name <- 1
bsubset_df <- btracks[btracks$seqnames == 1,]
seq_start <- which(abs(bsubset_df$start - 29000527) == min(abs(bsubset_df$start - 29000527)))
seq_end <- which(abs(bsubset_df$start - 29008920) == min(abs(bsubset_df$start - 29008920)))
bsubset_df <- bsubset_df[seq_start:seq_end,]
subset_df$exp <- 'clf'
asubset_df$exp <- 'clf-alp2'
bsubset_df$exp <- 'clf-alp1'
doub <- rbind(subset_df, asubset_df,bsubset_df)
p <- ggplot(doub, aes(x=start, y=score, fill=exp)) + 
    geom_bar(stat="identity", color="black", position=position_dodge(), size=0.1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


clf28 <- plot_peaks(bigwig = '~/clf28.bw',chrom = 4,start = 10382856-1000,end = 10388539+1000, genename = 'clf')
clf28_alp1_1 <- plot_peaks(bigwig = '~/clf28_alp1-1.bw',chrom = 4,start = 10382856-1000,end = 10388539+1000,genename = 'clf alp1-1')
clf28_alp2_1 <- plot_peaks(bigwig = '~/clf28_alp2-1.bw',chrom = 4,start = 10382856-1000,end = 10388539+1000,genename = 'clf alp2-1')
p <- ggpubr::ggarrange(clf28,clf28_alp1_1,clf28_alp2_1,nrow = 3,ncol =1 )
ggsave('~/Desktop/testplot.png', p)


clf28_alp1_chip_up <- dplyr::filter(clf28_alp1_chip, log2FoldChange > 1.5)
clf28_alp2_chip_up <- dplyr::filter(clf28_alp2_chip, log2FoldChange > 1.5)
clf28_alp1_rna_up <- dplyr::filter(clf28_alp1_rna, log2FoldChange > 1.5)
clf28_alp2_rna_up <- dplyr::filter(clf28_alp2_rna, log2FoldChange > 1.5)

clf28_alp1_chip_down <- dplyr::filter(clf28_alp1_chip, log2FoldChange < 0.75)
clf28_alp2_chip_down <- dplyr::filter(clf28_alp2_chip, log2FoldChange < 0.75)
clf28_alp1_rna_down <- dplyr::filter(clf28_alp1_rna, log2FoldChange < -1.5)
clf28_alp2_rna_down <- dplyr::filter(clf28_alp2_rna, log2FoldChange < -1.5)


up_genes <- list('clf28_alp1_chip_up'=clf28_alp1_chip_up$gene_id,'clf28_alp2_chip_up'=clf28_alp2_chip_up$gene_id,'clf28_alp1_rna_up'=clf28_alp1_rna_up$gene_id,'clf28_alp2_rna_up'=clf28_alp2_rna_up$gene_id)
chip_up_rna_down_genes <- list('clf28_alp1_chip_up'=clf28_alp1_chip_up$gene_id,'clf28_alp2_chip_up'=clf28_alp2_chip_up$gene_id,'clf28_alp1_rna_down'=clf28_alp1_rna_down$gene_id,'clf28_alp2_rna_down'=clf28_alp2_rna_down$gene_id)
chip_down_rna_up_genes <- list('clf28_alp1_chip_down'=clf28_alp1_chip_down$gene_id,'clf28_alp2_chip_down'=clf28_alp2_chip_down$gene_id,'clf28_alp1_rna_up'=clf28_alp1_rna_up$gene_id,'clf28_alp2_rna_up'=clf28_alp2_rna_up$gene_id)
down_genes <- list('clf28_alp1_chip_down'=clf28_alp1_chip_down$gene_id,'clf28_alp2_chip_down'=clf28_alp2_chip_down$gene_id,'clf28_alp1_rna_down'=clf28_alp1_rna_down$gene_id,'clf28_alp2_rna_down'=clf28_alp2_rna_down$gene_id)
overlaps_and_venn_output(genes=down_genes,out_file_prefix = '/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/clf28_as_control_RNA_CHIP_DOWN_overlaps',plot_title = 'clf28 as control RNA CHIP DOWN overlaps')


#Looking for motifs in the dataset where the IDs for ChIP and RNA-seq match
clf28_alp1_motif_df <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_transposons.csv') %>% dplyr::filter(gene_id %in% full$gene_id)
write.table(clf28_alp1_motif_df,file='/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/motif_analysis/clf28_alp1_rna_overlaps_motif.bed', sep = '\t', quote=F, row.names = F)
clf28_alp2_motif_df <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_transposons.csv') %>% dplyr::filter(gene_id %in% full$gene_id)
write.table(clf28_alp2_motif_df,file='/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/rna_chip_comparison/motif_analysis/clf28_alp2_rna_overlaps_motif.bed', sep = '\t', quote=F, row.names = F)






###clf28 double mutants overlap with single mutants
overlaps_double <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/differential_peaks/clf28_double_mutant_differential_peaks_sig_genes_overlaps.csv')
alp_specific_genes <- overlaps_double$clf28_alp1_diff_sig.clf28_alp2_diff_sig[grep('AT',overlaps_double$clf28_alp1_diff_sig.clf28_alp2_diff_sig)]
clf_control_alp1 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp1_me3_diff_v_clf28_annotated_both.csv')%>% filter(FC_KO > 1.5 | FC_KO > 1.5)
clf_control_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/2019_ChIP/epic2/clf_as_control/clf28_alp2_me3_diff_v_clf28_annotated_both.csv')%>% filter(FC_KO > 1.5 | FC_KO > 1.5)
genelist <- list(clf28_control_alp1=clf_control_alp1$feature,clf28_control_alp2=clf_control_alp2$feature,alp_specific_genes=alp_specific_genes)
overlaps_and_venn_output(genes=genelist,out_file_prefix='/Volumes/sesame/ALP_Omics/ChIP/validations/clf_double_mutants/clf_control/clf_alp_raw_v_clf_control_alp',plot_title='Genes specifc to alp double mutants in clf28 and clf28 double mutants relative to Col-0 vs. clf28 double mutants relative to clf28')
