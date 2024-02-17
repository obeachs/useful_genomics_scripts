library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gplots)
library(RVenn)
library(ggrepel)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(tidyverse)
library(coriell)
library(scales)
library(gtools)
library(writexl)
scale_fill_manual(values = c("H3K27me3 Down" = "
", "H3K27me3 Up" = "#40B0A6"))

theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)

salba_palette <- c('BUD Stage 6-7'='#ddd2f4','BUD Stage 8-9'='#b696f7',
 'BUD Stage 10-11'='#9463f7',
'GYN Stage 8-9'="#A6CEE3",'GYN Stage 10-11'="#1F78B4", 
'STA Stage 8-9'="#FDBF6F", 'STA Stage 10-11'="#FF7F00",
'PET Stage 8-9'="#faf8a4",'PET Stage 10-11'="#ECC846",
'SEP Stage 8-9'="#B2DF8A",'SEP Stage 10-11'="#33A02C")
names(salba_palette)
as.vector(salba_palette)
salba_palette_df <- data.frame('sample'=names(salba_palette), 'code'=as.vector(salba_palette))


plot_fold_changes <- function(gene, outprefix, se_done='TRUE'){
    x_axis_order <- c(
    "GYN Stage 8-9", "GYN Stage 10-11",
    "STA Stage 8-9", "STA Stage 10-11",
    "PET Stage 8-9", "PET Stage 10-11",
    "SEP Stage 8-9", "SEP Stage 10-11"
    )
    input_df <- read.csv('~/Salba_rna/genelists//organ_FPKMS.csv')
    organ_list <- list('GYN','STA','PET','SEP')
    df <- input_df %>% filter(grepl(gene, gene_id)) %>% 
        filter(!grepl('BUD', sample)) %>% 
            mutate(sample=gsub('_',' ', sample)) %>% mutate(sample=gsub('Stage 10', 'Stage 10-11', sample)) %>%
            mutate(sample=gsub('Stage 8', 'Stage 8-9', sample)) %>% mutate(sample=gsub('Stage 6', 'Stage 6-7', sample))
    print(df)
    # Create a bar plot with error bars
    p <- ggplot(df, aes(x = sample, y = avg_FPKM, fill = sample)) +
        geom_bar(stat = "identity", color='black',position = "dodge") +
        geom_errorbar(aes(ymin = avg_FPKM - (sqrt(sd_FPKM)/4), ymax = avg_FPKM + (sqrt(sd_FPKM))/4), position = position_dodge(0.9), width = 0.25) +
        labs(title = paste0('Average FPKM ', gene),
            x = "Sample",
            y = "Average FPKM")+
        guides(fill=guide_legend(title="Sample"))+
          theme(text = element_text(size = 20))+
        theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2))+
        theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"),
        panel.background = element_blank())+
        scale_fill_manual(values = salba_palette,limits = x_axis_order,
                    breaks = x_axis_order,
                    labels = x_axis_order)+
        scale_x_discrete(limits = x_axis_order)+
        theme

    ggsave(paste0(outprefix,gene,'_expression_fpkm.pdf'), height=10,width = 10, dpi=500)
    ggsave(paste0(outprefix,gene,'_expression_fpkm.png'), height=10,width = 10, dpi=500)
    }
plot_fold_changes('Sal02g14170L', outprefix = '~/thesis_figs_and_tables/salba/qpcr_rna-seq_comparisons/')


qpcr_plots <- function(anaylsis_file, out_prefix, organs=TRUE){
    data <- read.table(anaylsis_file, sep='\t', header=T) %>%
        mutate(sam_type= gsub('EARLY',' Stage 6-7',sam_type)) %>%
        mutate(sam_type= gsub('MID',' Stage 8-9',sam_type)) %>%
        mutate(sam_type= gsub('LATE',' Stage 10-11',sam_type))
    genes <- unique(data$sam_name)
    print(genes)
    x_axis_order <- c(
    "GYN Stage 8-9", "GYN Stage 10-11",
    "STA Stage 8-9", "STA Stage 10-11",
    "PET Stage 8-9", "PET Stage 10-11",
    "SEP Stage 8-9", "SEP Stage 10-11"
    )
    count <- 0
    seen <- list()
    for(gene in genes){
        count <- count + 1
        print(gene)
        gene <- gsub(' ','', gene)
        slim_data <- data %>% filter(grepl(gene,sam_name ))
        slim_data <- slim_data %>% mutate(cp=ifelse(cp > 100, log(cp), cp))%>% 
        mutate(sd=ifelse(sd > 100, log(sd), sd)) %>% filter(grepl(gene,sam_name )) %>% mutate(cp=ifelse(cp > 1000, (log(cp)/10)/10, cp))%>% 
        mutate(sd=ifelse(sd > 1000, (log(sd)/10), sd))
        print(slim_data)
        if(organs==TRUE){
        p <- ggplot(slim_data, aes(x=sam_type, y=cp, fill=sam_type)) + 
         geom_bar(stat="identity", color="black", position=position_dodge()) +
        geom_errorbar(aes(ymin=cp-(0.5)*sd, ymax=cp+(0.5)*sd), width=.2,position=position_dodge(.9))+
        theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"),
        panel.background = element_blank())+
        scale_fill_manual(values = salba_palette,limits = x_axis_order,
                    breaks = x_axis_order,
                    labels = x_axis_order)+
        scale_x_discrete(limits = x_axis_order)+
        labs(title = paste0('RT-qPCR ', gene),
            x = "Sample",
            y = "Value relative to reference TUB1")+
        guides(fill=guide_legend(title="Sample"))+
        theme(text = element_text(size = 20))+
        theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2))+theme+
          scale_x_discrete(limits = c(
        "GYN Stage 8-9", "GYN Stage 10-11",
        "STA Stage 8-9", "STA Stage 10-11",
        "PET Stage 8-9", "PET Stage 10-11",
        "SEP Stage 8-9", "SEP Stage 10-11"
        ))
        #      guides(fill = guide_legend(override.aes = list(order = c(
        # "GYN Stage 8-9", "GYN Stage 10-11",
        # "STA Stage 8-9", "STA Stage 10-11",
        # "PET Stage 8-9", "PET Stage 10-11",
        # "SEP Stage 8-9", "SEP Stage 10-11"
        # ))))
        ggsave(paste0(out_prefix,gene,'_',count,'_qpcr.pdf'),p,height = 10, width = 10, dpi=500)
    }
    }
}
qpcr_plots('~/Documents/Salba_qPCRs/JB_BRAPA_ORGANS_MID_LATE_307-311-317_qPCR_analysis.tsv','~/thesis_figs_and_tables/salba/qpcr_rna-seq_comparisons/' )

