library(CSAR)
library(org.At.tair.db)
library(GenomeGraphs)
library(bPeaks)
library(ggplot2)
library(ggpubr)
library(plotgardener)
library(RColorBrewer)
theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(12)

plot_peaks <- function(bigwig_list,chrom,start,end,genename,title_list, out){
     set_y <- 0
     df_list <- list()
     plot_list <- list()

     for(i in bigwig_list){
          print(i)
          tracks <- as.data.frame(plotgardener::readBigwig(i))
          tracks$seqnames <- gsub('Chr','', tracks$seqnames)
          chrom_name <- chrom
          subset_df <- tracks[tracks$seqnames == chrom,]
          seq_start <- which(abs(subset_df$start - start) == min(abs(subset_df$start - start)))
          seq_end <- which(abs(subset_df$start - end) == min(abs(subset_df$start - end)))
          subset_df <- subset_df[seq_start:seq_end,]
          subset_df$score[subset_df$score < 0] <- (subset_df$score[subset_df$score < 0])/-10
          df_list[[i]] <- subset_df
          max_h <- max(subset_df$score)
          if(max_h>set_y)
               {set_y <- max_h}
          print(max_h)
          print(set_y)
     }
     for(i in (1:length(df_list))){
          if(!grepl('Col', title_list[i])){
          title <- paste('Chr',chrom,sep = '')
          title <- paste(genename, title, sep = ' ')
          item <- as.data.frame(df_list[i])
          colnames(item) <- c('seqnames','start','end','width','strand','score')
          plot <- ggplot(data =item, aes(x= start ,y = score))+geom_area(fill =  mycolors[i]) +
          theme(text = element_text(size = 20))+
          scale_x_continuous(name = title, breaks = scales::pretty_breaks(n = 5))+
          scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+
          scale_y_continuous(name = 'Coverage/Peak Score')+
          scale_colour_manual(values = mycolors[i])+
          theme+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          ggtitle(title_list[i])+
          theme(plot.title = element_text(face = "italic"))+
          ylim(0, set_y)
          plot_list[[i]] <- plot
          }
          else{
          title <- paste('Chr',chrom,sep = '')
          title <- paste(genename, title, sep = ' ')
          item <- as.data.frame(df_list[i])
          colnames(item) <- c('seqnames','start','end','width','strand','score')
          plot <- ggplot(data =item, aes(x= start ,y = score))+geom_area(fill =  mycolors[9]) +
          theme+
          theme(text = element_text(size = 20))+
          scale_x_continuous(name = title, breaks = scales::pretty_breaks(n = 5))+
          scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+
          scale_y_continuous(name = 'Coverage/Peak Score')+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          ggtitle(title_list[i])+
          ylim(0, set_y)
          plot_list[[i]] <- plot
          }
     }
     ggsave(out,ggpubr::ggarrange(plotlist=plot_list, ncol = 1),width = 20, height = 15, limitsize = FALSE)
}



bigwig_clf = '/Volumes/sesame/joerecovery/Project_folder/alp_omics/Bennett_2019/bigwigs/CHKPEI85219060189-4_190716_X602_FCH2TJYCCX2_L4_CHKPEI85219060189-4_phix_removed_sorted.bw'
bigwig_clf_alp1 = '/Volumes/sesame/joerecovery/Project_folder/alp_omics/Bennett_2019/bigwigs/CHKPEI85219060188-10_190716_X602_FCH2TJYCCX2_L3_CHKPEI85219060188-10_phix_removed_sorted.bw'
bigwig_clf_alp2 = '/Volumes/sesame/joerecovery/Project_folder/alp_omics/Bennett_2019/bigwigs/CHKPEI85219060188-4_190716_X602_FCH2TJYCCX2_L3_CHKPEI85219060188-4_phix_removed_sorted.bw'
my_list <- list(bigwig_clf,bigwig_clf_alp1,bigwig_clf_alp2)
titles <- list('clf28','clf28 alp1-1', 'clf28 alp2-1')



alp2 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/alp2_rpkm.bw'
ref6 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/20cd21_ChIP_bams/ref6_rpkm.bw'
ref6_alp2 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/ref6_alp2_rpkm.bw'
col <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/col_rpkm.bw'
glist <- list(alp2, ref6, ref6_alp2, col)
titles <- list('alp2', 'ref6', 'ref6 alp2','Col-0')
true_ref6_alp2 <- read.csv('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/true_peaks/true_peaks.csv')
true_ref6_alp2_UP <- true_ref6_alp2[true_ref6_alp2$FC_KO>1,] %>% dplyr::distinct()
true_ref6_alp2_DOWN <- true_ref6_alp2[true_ref6_alp2$FC_KO<1,] %>% dplyr::distinct()
for(i in 1:nrow(true_ref6_alp2_DOWN)){
     chr <- as.numeric(gsub('Chr','',true_ref6_alp2_DOWN$chr[i]))
     start <- as.numeric(gsub('Chr','',true_ref6_alp2_DOWN$start[i]))
     end <- as.numeric(gsub('Chr','',true_ref6_alp2_DOWN$end[i]))
     feature <- true_ref6_alp2_DOWN$feature[i]
     outname <- paste0('/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_control/true_peaks/peaks_H3K27me3_DOWN/', feature,'.png')
     plot_peaks(glist, chr, start-2000, end+2000, feature, titles, outname)
}
plot_peaks(glist,3,2768688-2000,2771174 + 2000,'AT3G09070', titles, '/Volumes/sesame/ALP_Omics/ChIP/validations/ref6_alp2/ref6_top_peaks/AT3G09070.png')

cralp2 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2021_ChIP_bams/alp2-1_cpm_bamcompare.bw'
bt_alp2 <- '/Volumes/sesame/joerecovery/Project_folder/alp_omics/2020_ChIP_bams/alp2-1_cpm_bamcompare.bw'
glist <- list(cralp2,bt_alp2)
titles <- list('alp2_2020','alp2_2021')
plot_peaks(glist,5,3173080-5000,3179752+5000,'FLC', titles)


plot_peaks(list(bigwig_clf,bigwig_clf_alp1,bigwig_clf_alp2),5,3173080,3179752,'FLC', titles)


clf <- plot_peaks(bigwig_clf,4,10382572-20000,10388823 + 20000,'AG')
clf_alp1 <- plot_peaks(bigwig_clf_alp1,4,10382572-20000,10388823 + 20000,'AG')
clf_alp2 <- plot_peaks(bigwig_clf_alp2,4,10382572-20000,10388823 + 20000,'AG')
plot_list <- list(clf,clf_alp1,clf_alp2)
ggpubr::ggarrange(clf,clf_alp1,clf_alp2, ncol = 1)


tracks <- as.data.frame(plotgardener::readBigwig(bigwig_cpm))
subset_df <- tracks[tracks$seqnames == 'Chr4',]
seq_start <- which(abs(subset_df$start - (10382572-20000)) == min(abs(subset_df$start - (10382572-20000))))
seq_end <- which(abs(subset_df$start - (10388823+20000)) == min(abs(subset_df$start - (10388823+20000))))
subset_df <- subset_df[seq_start:seq_end,]

