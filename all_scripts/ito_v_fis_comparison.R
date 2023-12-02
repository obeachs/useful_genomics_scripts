source('/Volumes/sesame/joerecovery/scripts/formattable_functions_and_tables.R')
source('~/Salba_RNA/scripts/GO_dotplotter.R')
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(RColorBrewer)
library(webr)
library(cowplot)
library(ggrepel)
library(coriell)
  theme <- theme(
    axis.text.x = element_text(colour = "black"),
    panel.background = element_blank(), panel.border = element_rect(fill = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "line")
  )

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

panther_go_maker <- function(genelist,out_name){
  GO_out <- panther_go(gene_list = genelist,organism = 3702,annot_dataset = 'biological_process')
  print(head(GO_out))
  GG <- GO_out %>% dplyr::select(term,fold_enrichment, fdr) %>% filter(fdr < 0.05)
  GG <- GG[!grepl('GO:', GG$term),] %>% arrange(fdr)%>% mutate(term=unlist(term))
  write.csv(GG,out_name, quote=F, row.names=F)
  }

symbols <- read.csv('~/Salba_RNA/genelists/all_hits_names.csv') %>% dplyr::select(tair, gene.name)
ito <- read.table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/dex_v_mock_35S-SUP-GR_DEGs.tsv', header=T)%>% filter(adj.P.Val <0.05) %>% 
dplyr::rename(tair=ORF) %>% left_join(symbols) %>% dplyr::select(tair, gene.name, logFC) %>% dplyr::rename('Gene ID'=tair, Symbol=gene.name,  log2FC=logFC) %>% distinct()
make_nice_table(ito, colnums=c(3), '~/thesis_figs_and_tables/sup/ito_significant_DEGs.png')


uppydowny<-data.frame( Status=c('Up', 'Down'),Count=c(length(ito$ID[ito$logFC>1]),length(ito$ID[ito$logFC<1])))
uppydowny$ratio <- round((uppydowny$Count/sum(uppydowny$Count))*100)
ito_for_GO <- ito %>% filter(adj.P.Val < 0.05)
panther_go_maker(ito_for_GO$ORF, '~/thesis_figs_and_tables/sup/ito_qval_0-05_pantherGO')
go_dotplotter('~/thesis_figs_and_tables/sup/ito_qval_0-05_pantherGO')
panther_go_maker(ito$ORF, '~/thesis_figs_and_tables/sup/ito_pval_0-05_pantherGO')
go_dotplotter('~/thesis_figs_and_tables/sup/ito_pval_0-05_pantherGO')


raw <- ggplot(uppydowny, aes(x = 2, y = Count, fill = Status)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = c("#f54d4d","#a8d671")) +
  theme_void()+
  xlim(0.5, 2.5)+
  annotate(geom = 'text', x = 0.5, y = 0,size=10, label = paste0('Total: ',sum(uppydowny$Count)))+
   theme(legend.text = element_text(size=30))+
   theme(legend.title = element_text(size=30))+
    geom_text(aes(label = paste0(ratio,'%')),
            position = position_stack(vjust = 0.5),size=10, color = "white") +
  theme(axis.text = element_text(size = 10, colour="white"))










make_nice_table(uppydowny, outname='~/thesis_figs_and_tables/sup/ito_up_v_DOWN_q0-05.png')
write.csv(uppydowny,'~/thesis_figs_and_tables/sup/ito_up_v_DOWN_q0-05.csv', row.names=F, quote=F)
ggsave('~/thesis_figs_and_tables/sup/ito_up_v_DOWN_q0-05.pdf', raw, height = 10, width=10)

mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(2)
fis_2 <- read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/merge2.5D_sorted_0.05.txt') %>% dplyr::select(locus, Symbols, logFC) %>%
dplyr::rename(Symbol=Symbols, log2FC=logFC) %>% mutate(Symbol=ifelse(is.na(Symbol),locus, Symbol )) %>% mutate(Symbol=ifelse(Symbol=='',locus, Symbol )) %>% dplyr::rename('Gene ID'=locus) %>% distinct() %>% 
arrange(desc(abs(log2FC))) %>% na.omit()
fis_3 <-read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/merge3.5D_sorted_0.05.txt') %>% dplyr::select(locus, Symbols, logFC) %>%
dplyr::rename(Symbol=Symbols, log2FC=logFC) %>% mutate(Symbol=ifelse(is.na(Symbol),locus, Symbol )) %>% mutate(Symbol=ifelse(Symbol=='',locus, Symbol )) %>% dplyr::rename('Gene ID'=locus) %>% distinct() %>%
arrange(desc(abs(log2FC))) %>% slice_head(n=30)
fis_5 <-read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/merge5D_sorted_0.05.txt') %>% dplyr::select(locus, Symbols, logFC) %>%
dplyr::rename(Symbol=Symbols, log2FC=logFC) %>% mutate(Symbol=ifelse(is.na(Symbol),locus, Symbol )) %>% mutate(Symbol=ifelse(Symbol=='',locus, Symbol )) %>% dplyr::rename('Gene ID'=locus) %>% distinct() %>%
arrange(desc(abs(log2FC)))%>% slice_head(n=30)


make_nice_table(fis_5, colnums = c(3), outname = '~/thesis_figs_and_tables/sup/5DAI_genes.png')


overlaps_and_venn_output(genes=list('2.5 DAI'=fis_2$locus, '35S:SUP-GR'=ito$ORF), out_file_prefix = '~/thesis_figs_and_tables/sup/ito_fis2.5_overlaps',
 'Overlap of Significant DEGs 35S::SUP-GR and FIS sup-5 2.5 DAI')
 overlaps_and_venn_output(genes=list('3.5 DAI'=fis_3$locus, '35S:SUP-GR'=ito$ORF), out_file_prefix = '~/thesis_figs_and_tables/sup/ito_fis3.5_overlaps',
 'Overlap of Significant DEGs 35S::SUP-GR and FIS sup-5 3.5 DAI')
overlaps_and_venn_output(genes=list('5 DAI'=fis_5$locus, '35S:SUP-GR'=ito$ORF), out_file_prefix = '~/thesis_figs_and_tables/sup/ito_fis5_overlaps',
 'Overlap of Significant DEGs 35S::SUP-GR and FIS sup-5 5 DAI')

itofis2 <- ggVennDiagram::ggVennDiagram(list('2.5 DAI'=fis_2$locus, '35S:SUP-GR'=ito$ORF)) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='FIS sup-5 2.5 DAI') +  scale_color_brewer(palette = "Paired")+
    theme(text = element_text(size = 50),plot.title = element_text(size = 16),legend.title = element_text(size = 12))+
    guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+theme(plot.margin = margin(0, 0, 0, 0))
itofis3 <- ggVennDiagram::ggVennDiagram(list('3.5 DAI'=fis_3$locus, '35S:SUP-GR'=ito$ORF)) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='FIS sup-5 3.5 DAI') +  scale_color_brewer(palette = "Paired")+
    theme(
    text = element_text(size = 50),          # Adjust the font size for text
    plot.title = element_text(size = 16),    # Adjust the title font size
    legend.title = element_text(size = 12)  # Adjust the legend title font size
  )+
   guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+
   theme(plot.margin = margin(0, 0, 0, 0))
itofis5 <- ggVennDiagram::ggVennDiagram(list('5 DAI'=fis_5$locus, '35S:SUP-GR'=ito$ORF)) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='FIS sup-5 5 DAI') +  scale_color_brewer(palette = "Paired")+
    theme(
    text = element_text(size = 50),          # Adjust the font size for text
    plot.title = element_text(size = 16),    # Adjust the title font size
    legend.title = element_text(size = 12)  # Adjust the legend title font size
  )+
   guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+
   theme(plot.margin = margin(0, 0, 0, 0))



b <- plot_grid(
  plot_grid(itofis2, itofis3, nrow = 1, ncol = 2),
  plot_grid(NULL, itofis5, NULL, nrow = 1, rel_widths = c(0.5, 1, 0.5)),
  nrow = 2
)
ggsave('~/thesis_figs_and_tables/sup/ito_qvalue_0.05_fis_qvalue_0.05_overlaps.pdf', b, height = 10, width = 10, dpi = 4000)
















ito_2 <- read.table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/dex_v_mock_35S-SUP-GR_DEGs.tsv',header=T) %>% 
filter(P.Value <0.05) %>% dplyr::rename('locus'=ORF) %>% inner_join(fis_2, by='locus') %>% dplyr::select(locus,logFC.x,logFC.y)
ito_3 <- read.table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/dex_v_mock_35S-SUP-GR_DEGs.tsv',
header=T) %>% filter(P.Value <0.05) %>% dplyr::rename('locus'=ORF) %>% inner_join(fis_3, by='locus')%>% dplyr::select(locus,logFC.x,logFC.y)
ito_5 <- read.table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/dex_v_mock_35S-SUP-GR_DEGs.tsv',
header=T) %>% filter(P.Value <0.05) %>% dplyr::rename('locus'=ORF) %>% inner_join(fis_5, by='locus')%>% dplyr::select(locus,logFC.x,logFC.y)


reshaped_data <- ito_5 %>%
  tidyr::gather(logFC_type, logFC_value, logFC.x, logFC.y) %>% mutate(logFC_type=ifelse(logFC_type=='logFC.x', '35S:SUP-GR','FIS sup-5'))
result <- reshaped_data %>%
  group_by(locus, logFC_type) %>%
  slice(which.max(abs(logFC_value))) %>%
  ungroup()
p<-ggplot(result, aes(x = locus, y = logFC_value, fill=logFC_type)) +
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  #geom_errorbar(width = 0.2, position = position_dodge(0.9)) +
  labs(
       x = "TAIR ID",
       y = "log2FC")+
       theme+
  # scale_fill_brewer(palette = mycolors)+  
  scale_fill_manual(values = mycolors)+
  ylim(min(reshaped_data$logFC_value) -1, max(reshaped_data$logFC_value) +1)+
  geom_hline(yintercept = 0)+
  # facet_wrap(~face)+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"))+
  theme(text = element_text(size = 20))+
  #theme(legend.position="none") +
  #scale_y_continuous(breaks = pretty(df_check$TPM, n = 10))+
  scale_y_continuous(n.breaks=20)+
  guides(fill=guide_legend(title="Tissue and Stage"))+
  theme+
  theme(axis.title.x = element_text(margin = margin(t = 20)),
          axis.title.y = element_text(margin = margin(r = 20)),
          axis.text.x = element_text(margin = margin(t = 10)),
          axis.text.y = element_text(margin = margin(r = 10)))+
  ggtitle('Differential expression of shared genes - 35S::SUP-GR and FIS sup-5 microarrays')

ggsave('~/thesis_figs_and_tables/sup/ito_fis5_expresion_barplot.pdf', p, width=10, height=11)



up2 <- fis_2[fis_2$logFC>0,]
down2 <- fis_2[fis_2$logFC<0,]

up3 <- fis_3[fis_3$logFC>0,]
down3 <- fis_3[fis_3$logFC<0,]

up5 <- fis_5[fis_5$logFC>0,]
down5 <- fis_5[fis_5$logFC<0,]


head(ito)

allgeneslist <- list(fis_2$locus,fis_3$locus,fis_5$locus)
downgeneslist <- list(down2$locus, down3$locus, down5$locus)
upgeneslist <- list(up2$locus, up3$locus, up5$locus)




library(cowplot)
library(ggplot2)


names(allgeneslist) <- c('2.5 DAI', '3.5 DAI', '5 DAI')
names(upgeneslist) <- c('2.5 DAI', '3.5 DAI', '5 DAI')
names(downgeneslist) <- c('2.5 DAI', '3.5 DAI', '5 DAI')
overlaps_and_venn_output(genes=allgeneslist, out_file_prefix='~/thesis_figs_and_tables/sup/microarray_DEGS_comparison', plot_title='Overlap of all DEGs at 2.5, 3.5 and 5 DAI')
overlaps_and_venn_output(genes=downgeneslist, out_file_prefix='~/thesis_figs_and_tables/sup/microarray_down_DEGS_comparison', plot_title='Overlap of downregulated DEGs at 2.5, 3.5 and 5 DAI')
overlaps_and_venn_output(genes=upgeneslist, out_file_prefix='~/thesis_figs_and_tables/sup/microarray_up_DEGS_comparison', plot_title='Overlap of upregulated DEGs at 2.5, 3.5 and 5 DAI')
allp <- ggVennDiagram::ggVennDiagram(allgeneslist) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='All Significant DEGs') +  scale_color_brewer(palette = "Paired")+
    theme(text = element_text(size = 50),plot.title = element_text(size = 16),legend.title = element_text(size = 12))+
    guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+theme(plot.margin = margin(0, 0, 0, 0))
downp <- ggVennDiagram::ggVennDiagram(downgeneslist) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='Significant downregulated DEGs') +  scale_color_brewer(palette = "Paired")+
    theme(
    text = element_text(size = 50),          # Adjust the font size for text
    plot.title = element_text(size = 16),    # Adjust the title font size
    legend.title = element_text(size = 12)  # Adjust the legend title font size
  )+
   guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+
   theme(plot.margin = margin(0, 0, 0, 0))
upp <- ggVennDiagram::ggVennDiagram(upgeneslist) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkorchid",)+
    ggplot2::scale_x_continuous(expand = expansion(mult = .2))+
    labs(title ='Signficiant upregulated DEGs') +  scale_color_brewer(palette = "Paired")+
    theme(
    text = element_text(size = 50),          # Adjust the font size for text
    plot.title = element_text(size = 16),    # Adjust the title font size
    legend.title = element_text(size = 12)  # Adjust the legend title font size
  )+
   guides(fill = 'none')+theme(plot.title = element_text(hjust = 0.5))+
   theme(plot.margin = margin(0, 0, 0, 0))


ggpubr::ggarrange(allp, downp,upp, nrow=1)
b <- plot_grid(
  plot_grid(upp, downp, nrow = 1, ncol = 2),
  plot_grid(NULL, allp, NULL, nrow = 1, rel_widths = c(0.5, 1, 0.5)),
  nrow = 2
)
ggsave('~/thesis_figs_and_tables/sup/microarray_DEGS_comparison_ggpubr.pdf', b, height = 10, width = 10,dpi = 400000)
# Arrange the top two plots
top_arrangement <- ggpubr::ggarrange(
  upp + theme_void(),
  downp + theme_void(),
  ncol = 2
)
# Create a blank plot for the bottom plot area
blank_plot <- ggplot() + theme_void()
# Arrange the top two plots and the blank plot
final_arrangement <- cowplot::plot_grid(
  top_arrangement,
  blank_plot,
  allp + theme_void(),
  ncol = 1,
  axis = "tblr"
)
# Plot the final arrangement
ggsave('~/thesis_figs_and_tables/sup/microarray_DEGS_comparison_ggpubr.pdf',final_arrangement, width = 20, height=15)
















DAI <- c('Chr','Chr2','Chr3','Chr4','Chr5')
count <- c(lchr1,lchr2,lchr3,lchr4,lchr5)
hicup_count <- c(3782,2440,2928,2318,3294)
df <- tibble(DAI=c('0 DAI','2.5 DAI', '3.5 DAI','5 DAI'),
Downregulated=c(0,length(down2$locus), length(down3$locus), length(down5$locus)),
Upgregulated=c(0,length(up2$locus), length(up3$locus), length(up5$locus)),
)
write.csv(df, '~/thesis_figs_and_tables/sup/microarray_DEG_counts_table2.csv', row.names=F, quote=F)

make_nice_table(df,outname='~/thesis_figs_and_tables/sup/microarray_DEG_counts_table2.png')

samps <- c('2.5 DAI','2.5 DAI','3.5 DAI','3.5 DAI','5 DAI','5 DAI')
condition <- c('Up','Down','Up','Down','Up','Down')
value <- c(length(up2$locus),length(down2$locus),length(up3$locus),length(down3$locus),length(up5$locus),length(down5$locus))
data <- data.frame(samps,condition,value)
data <- read.csv('~/thesis_figs_and_tables/sup/microarray_DEG_counts_table1.csv')

 p <- ggplot(data, aes(fill=condition, y=value, x=samps,label = value)) + 
      theme +
      geom_bar(position="stack", stat="identity") +
      geom_text(size = 8, position = position_stack(vjust = 0.5),colour = "white")+
      scale_fill_manual(values=c("#f54d4d", "#a8d671"))+
      scale_x_discrete(limits = c("2.5 DAI", "3.5 DAI", "5 DAI"))+
      ggtitle("Differentially expressed genes FIS sup-5 -/-  vs FIS sup-5 +/+") +
      guides(fill=guide_legend(title="DEG Status"))+
      theme(text = element_text(size = 20))+
      xlab("Days after DEX induction")+
      ylab("Gene Number")

ggsave('~/thesis_figs_and_tables/sup/microarray_DEG_counts.pdf',p, height=15, width=11)
