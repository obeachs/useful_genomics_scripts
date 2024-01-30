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
library(kableExtra)
library(tidyverse)
library(gggenomes)
library(wesanderson)
source('~/useful_genomics_scripts/all_scripts/formattable_functions_and_tables.R')


scale_fill_manual(values = c("H3K27me3 Down" = "#E1BE6A", "H3K27me3 Up" = "#40B0A6"))

theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)

wes_palette <- c(
'#DD8D29',
'#E2D200',
'#46ACC8',
'#E58601',
'#B40F20'
)
##QUAST results 
jgi <- read.csv('~/thesis_figs_and_tables/salba/JGI_quast.csv',header=T)
yang <- read.delim('~/thesis_figs_and_tables/salba/yang_quast.tsv', sep = '\t', header = T)
nanopore <- read.delim('~/thesis_figs_and_tables/salba/nanopore_SPAdes_assembly_quast.tsv', sep='\t', header = T)
tab <- rbind(jgi, yang, nanopore) %>% dplyr::select(X..contigs,Largest.contig,Total.length, 'GC....',N50, N90, L50,L90)
names(tab) <- c('Contigs','Largest Contig', 'Total length', 'GC Content', 'N50', 'N90','L50','L90') 
tab[] <- sapply(tab, as.numeric) 
tab <- as.data.frame(t(tab))
tab <- tab %>% mutate_all(~gsub("\\.00", "", as.character(.)))
colnames(tab)<- c('JGI Assembly', 'Yang et al 2023 Assembly', 'Nanopore reads assembly')

t <- kbl(tab, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down","add_linespace=-1.5mm")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:3, border_left = FALSE, border_right = FALSE)
cat(t)

# Arabidopssi stages
smyth <- read.csv('~/Salba_RNA/genelists/arabidopsis_stages.txt')
t <- kbl(smyth, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down","add_linespace=-1.5mm")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:2, border_left = FALSE, border_right = FALSE)
cat(t)


##Nanoplot stats
run_info <- read.csv('~/thesis_figs_and_tables/salba/nanopore_run_stats/run_stats', header=F)
t <- kbl(run_info, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down","add_linespace=-1.5mm")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:3, border_left = FALSE, border_right = FALSE)
cat(t)

### qPCR reference setup ----
  df <- data.frame('Tissue'=c('Bud', 'Whole inflorescence','Leaf','Bud', 'Whole inflorescence','Leaf'),
  'Mean Cp value'=c(22.5,22.1,22.4,17.6,17.9,18.1),'se'=c(0.9,0.7,0.65,1,1.2,1.21),'gene'=c('TUB1','TUB1','TUB1','GADPH','GADPH','GADPH'))
  p <- ggplot(df, aes(x=Tissue, y=Mean.Cp.value, fill=gene)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=Mean.Cp.value-(0.5)*se, ymax=Mean.Cp.value+(0.5)*se), width=.2,position=position_dodge(.9))+
    theme(axis.ticks.x = element_blank(),
    panel.background = element_blank())+
    guides(fill=guide_legend(title="Sample type"))+
    scale_fill_manual(values=c('TUB1'='#bedb8c','GADPH'='#f3f580'))+
    ylab('Mean cp value')+
    xlab('Gene name and Sample')+
    theme(text = element_text(size = 20))+
    theme(axis.text.x=element_text(hjust=0.95,vjust=0.2))+theme
    #axis.text.x = element_text(angle = 90)

    ggsave('~/thesis_figs_and_tables/salba/reference_qpcrs.pdf',p,height = 10, width = 10)

### transcript counts ----
  original <- read.table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/Sal.Chr.20210627.gff', 
  sep='\t', header=F) %>%  filter(V3=='gene') %>% filter(grepl('Contig',V1))
  original <- original %>% filter(V3=='gene') %>% filter(!grepl('Contig',V1))
  ref <- read.table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/merged.gtf', 
  sep='\t', header=F)
  ref <- ref %>% filter(V3=='transcript')
  gffcompare <- read.table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/compare.merged.gtf.tmap', sep='\t', header=T)
  u <- gffcompare %>% filter(class_code=='u') %>% filter(!grepl('Contig', qry_gene_id))
  transcript_data <- data.frame('Transcript classification'=c('Matches reference', 'Novel Gene', 'Novel transcript'), 'Count'=c(40623,3978, 4559))
  gffcompare_gtf <- read.table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/compare_opt.combined.gtf', sep='\t', header=F) %>% filter(grepl('Contig', V1))
  gffcompare_gtf <- gffcompare_gtf %>% filter(V3=='transcript') %>% filter(grepl('class_code u', V9))
  t <- kbl(transcript_data, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1),bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:ncol(transcript_data), border_left = FALSE, border_right = FALSE)
  cat(t)
##Alignment stats
stats <- read.csv('~/Salba_RNA/alignment_stats.csv') %>% mutate(Sample=gsub('_sorted','', Sample)) %>% mutate(Sample.1=gsub('_sorted','', Sample.1))
t <- kbl(stats, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1,3),bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:4, border_left = FALSE, border_right = FALSE)
cat(t)

##organ-specific genes
stage_8 <- data.frame(Organ=c('Gynoecium','Stamen','Petal','Sepal'),Count=c(1173,4099,469,1323))
t <- kbl(stage_8, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1),bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:4, border_left = FALSE, border_right = FALSE)
cat(t)
stage_10 <- data.frame(Organ=c('Gynoecium','Stamen','Petal','Sepal'),Count=c(522,3273,504,2742))
t <- kbl(stage_10, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1),bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:4, border_left = FALSE, border_right = FALSE)
cat(t)


# early and late flowering genes ---- 
  full_table <- read.csv('~/Salba_RNA/genelists/key_players.txt')
  early_genes <-  full_table %>% filter(time =='Early') %>% dplyr::select(gene.name)
  total <- nrow(early_genes)
  early_1 <- early_genes[1:(total/2),]
  early_2 <- early_genes[((total/2)+1):total,]
  early_new <- cbind(early_1,early_2)
  t <- kbl(early_new, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:ncol(early_new), border_left = FALSE, border_right = FALSE)

  late_genes <-  full_table %>% filter(time =='Late') %>% dplyr::select(gene.name)
  total <- nrow(late_genes)
  late_1 <- late_genes[1:(total/2),]
  late_2 <- late_genes[((total/2)+1):total,]
  late_new <-  cbind(late_1, late_2)
   t <- kbl(late_new, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:ncol(late_new), border_left = FALSE, border_right = FALSE)
##MADS homologs ---- 
species_counts <- read.csv('~/Salba_RNA/results/flower/mads/mads_homolog_counts.csv')
  species_counts <- read.csv('~/Salba_RNA/results/flower/mads/mads_homolog_counts.csv')
  names(species_counts) <- gsub('_',' ', names(species_counts))
  names(species_counts) <- gsub('\\.','-', names(species_counts))
  t <- kbl(species_counts, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
    column_spec(c(1,2),bold = TRUE) %>% 
    kable_styling() %>%
    column_spec(1:4, border_left = FALSE, border_right = FALSE)
  cat(t)


  tair_mads <- read.delim('~/Salba_RNA/genelists/MADS.csv', header = F, sep='\t')
  genes <- read.table('~/Salba_RNA/genelists/all_hits_strict.csv', header=T) %>% 
  filter(tair %in% tair_mads$V2)
  genes$gene_id <- gsub(",+$", "", genes$gene_id)
  t <- kbl(genes, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
    column_spec(c(1,2),bold = FALSE) %>% 
    kable_styling()
    column_spec(1:4, border_left = FALSE, border_right = FALSE)
  cat(t)
  nnew_gennes <- genes %>% separate_rows('gene_id')
  length(unique(nnew_gennes$gene_id))

#Trichome gene info 
trichome_grn_arabidopsis <- read.csv('~/Salba_RNA/genelists/trichome_GRN_salba.csv') %>% dplyr::select(gene_name, effect) %>%
  dplyr::rename(Gene_Name=gene_name,Effect=effect) %>% distinct()
  total_rows <- nrow(trichome_grn_arabidopsis)

  rows_per_part <- ceiling(total_rows / 3)

  df1 <- trichome_grn_arabidopsis[1:rows_per_part, ] %>% arrange(Effect)
  df2 <- trichome_grn_arabidopsis[(rows_per_part+1):(2*rows_per_part),]%>% arrange(Effect)
  df3 <- trichome_grn_arabidopsis[(2*rows_per_part + 1):total_rows, ]%>% arrange(Effect)
  df3[35,] <- df3[34,]
  df3[36,] <- df3[34,]
  combined_table <- cbind(df1, df2, df3)
  t <- kbl(combined_table, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1,3,5),bold = T) %>% 
  #kable_styling(latex_options = c("striped", "scale_down")) %>% 
  column_spec(1:ncol(combined_table), border_left = FALSE, border_right = FALSE)
  cat(t)







#Distribution of read lengths
jgi_tab <- read.table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/JOI_info/joi_readlength_counts.txt', header=T) %>%
 mutate(Data='PacBio JGI') %>% filter(Length > 0)
tab <- read.table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/nanopore_reads_raw/all_nanopores_combind_1000bp.txt', header=T) %>% 
mutate(Data='Nanopore Reads')%>% filter(Length > 0)
full <- rbind(jgi_tab,tab)


p <- ggplot(full, aes(x=Length, y=Count, fill=Data)) + 
    geom_bar(stat="identity", position=position_dodge(), size=0.1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme+
    scale_fill_manual(values=c("Nanopore Reads" = "#46ACC8", "PacBio JGI" = "#E2D200"))+
    theme(plot.title = element_text(hjust = 0.5),plot.margin = unit(c(1, 1, 1, 1), "line"))+
    xlab(" > Read length")+  
    ylab("Count")+
    scale_y_continuous(labels = scales::comma_format())+
    scale_x_continuous(expand = c(0, 0))

ggsave('~/Desktop/test_plot.pdf',p, height=10, width=10)
ggplot(tab, aes(x=Length, y=Count)) +
  geom_point() + 
  geom_segment( aes(x=Length, xend=Length, y=0, yend=Count))


##trichome GRN table
data <- read.csv('~/Salba_RNA/genelists/trichome_GRN_salba.csv') %>%
 mutate(tair=ifelse(tair!=gene_name,paste0(tair,'/',gene_name), tair))

new_data <- data %>% group_by(tair) %>%
  summarise(gene_ids = toString(unique(gene_id))) %>% arrange(tair)

total_rows <- nrow(new_data)

rows_per_part <- ceiling(total_rows / 3)

# Step 2: Split the table into three parts
df1 <- new_data[1:rows_per_part, ]
df2 <- new_data[(rows_per_part+1):(2*rows_per_part),]
df3 <- new_data[(2*rows_per_part + 1):total_rows, ]
combined_table <- bind_cols(df1, df2, df3)
t <- kbl(df1, format = "latex", booktabs = TRUE, latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(c(1),bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:ncol(combined_table), border_left = FALSE, border_right = FALSE)
cat(t)





# Differential expression analysis table 10-11 ---- 
    gyn_v_sta <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_STA_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 ) %>%
    gyn_v_pet <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_PET_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05 ) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
    gyn_v_sep <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_SEP_gtf_redo_no_lrt') %>%
      dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
    gyn_info <- data.frame('Comparison'=c('GYN vs STA','GYN vs PET', 'GYN v SEP'), 
    'Genes_Up'=c(nrow(gyn_v_sta %>% filter(log2FoldChange >0)),nrow(gyn_v_pet %>% filter(log2FoldChange >0)),nrow(gyn_v_sep %>% filter(log2FoldChange >0))),
    'Genes_Down'=c(-1*(nrow(gyn_v_sta %>% filter(log2FoldChange <0))),-1*(nrow(gyn_v_pet %>% filter(log2FoldChange <0))),-1*(nrow(gyn_v_sep %>% filter(log2FoldChange <0)))))
    names(gyn_info) <- gsub('\\.',' ', names(gyn_info))


    sta_v_gyn <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_STA_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2)
    sta_v_pet <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/STA_v_PET_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05 ) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
    sta_v_sep <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/sta_v_SEP_gtf_redo_no_lrt') %>%
      dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
    sta_info <- data.frame('Comparison'=c('STA vs GYN','STA vs PET', 'STA v SEP'), 
    'Genes_Up'=c(nrow(sta_v_gyn %>% filter(log2FoldChange < 0)),nrow(sta_v_pet %>% filter(log2FoldChange >0)),nrow(sta_v_sep %>% filter(log2FoldChange >0))),
    'Genes_Down'=c(-1*(nrow(sta_v_gyn %>% filter(log2FoldChange >0))),-1*(nrow(sta_v_pet %>% filter(log2FoldChange <0))),-1*(nrow(sta_v_sep %>% filter(log2FoldChange <0)))))
    names(gyn_info) <- gsub('\\.',' ', names(gyn_info))

    pet_v_gyn <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_PET_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2)
    pet_v_sta <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/STA_v_PET_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05 ) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
    pet_v_sep <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/PET_v_SEP_gtf_redo_no_lrt') %>%
      dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
    pet_info <- data.frame('Comparison'=c('PET vs GYN','PET vs STA', 'PET v SEP'), 
    'Genes_Up'=c(nrow(pet_v_gyn %>% filter(log2FoldChange < 0)),nrow(pet_v_sta %>% filter(log2FoldChange < 0)),nrow(pet_v_sep %>% filter(log2FoldChange >0))),
    'Genes_Down'=c(-1*(nrow(sta_v_gyn %>% filter(log2FoldChange >0))),-1*(nrow(pet_v_sta %>% filter(log2FoldChange >0 ))),-1*(nrow(pet_v_sep %>% filter(log2FoldChange <0)))))

    sep_v_gyn <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_SEP_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2)
    sep_v_sta <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/STA_v_SEP_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05 ) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
    sep_v_pet <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/PET_v_SEP_gtf_redo_no_lrt') %>%
      dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
    sep_info <- data.frame('Comparison'=c('SEP vs GYN','SEP vs STA', 'SEP v PET'), 
    'Genes_Up'=c(nrow(pet_v_gyn %>% filter(log2FoldChange < 0)),nrow(pet_v_sta %>% filter(log2FoldChange < 0)),nrow(pet_v_sep %>% filter(log2FoldChange <0))),
    'Genes_Down'=c(-1*(nrow(sta_v_gyn %>% filter(log2FoldChange >0))),-1*(nrow(pet_v_sta %>% filter(log2FoldChange > 0 ))),-1*(nrow(pet_v_sep %>% filter(log2FoldChange > 0)))))

  
  gyn_p <- ggplot(gyn_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    ylim(-5000, 5000)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Gynoecium')+
    ylab("Gene count")

  sta_p <- ggplot(sta_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    ylim(-5000, 5000)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Stamen')+
    ylab("Gene count")

    pet_p <- ggplot(pet_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    ylim(-5000, 5000)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Petal')+
    ylab("Gene count")  

  sep_p <- ggplot(gyn_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    ylim(-5000, 5000)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    )+ theme+
    ggtitle('Sepal')+
    ylab("Gene count")

  p <- ggpubr::ggarrange(gyn_p, sta_p, pet_p, sep_p, nrow = 1, ncol=4)
  ggsave('~/Salba_RNA/results/organs/tissue/stage_8-9/DEGs_lollipop.pdf', p,width = 210, height = 150, units = "mm", limitsize = F, dpi=500)


# DEGs Stage 8-9 -----
  gyn_v_sta <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_STA_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  gyn_v_pet <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_PET_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05 ) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  gyn_v_sep <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_SEP_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  gyn_info <- data.frame('Comparison'=c('GYN vs STA','GYN vs PET', 'GYN v SEP'), 
  'Genes_Up'=c(nrow(gyn_v_sta %>% filter(log2FoldChange >0)),nrow(gyn_v_pet %>% filter(log2FoldChange >0)),nrow(gyn_v_sep %>% filter(log2FoldChange >0))),
  'Genes_Down'=c(-1*(nrow(gyn_v_sta %>% filter(log2FoldChange <0))),-1*(nrow(gyn_v_pet %>% filter(log2FoldChange <0))),-1*(nrow(gyn_v_sep %>% filter(log2FoldChange <0)))))



  sta_v_gyn <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_STA_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2)
  sta_v_pet <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/STA_v_PET_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05 ) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  sta_v_sep <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/sta_v_SEP_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  sta_info <- data.frame('Comparison'=c('STA vs GYN','STA vs PET', 'STA v SEP'), 
  'Genes_Up'=c(nrow(sta_v_gyn %>% filter(log2FoldChange < 0)),nrow(sta_v_pet %>% filter(log2FoldChange >0)),nrow(sta_v_sep %>% filter(log2FoldChange >0))),
  'Genes_Down'=c(-1*(nrow(sta_v_gyn %>% filter(log2FoldChange >0))),-1*(nrow(sta_v_pet %>% filter(log2FoldChange <0))),-1*(nrow(sta_v_sep %>% filter(log2FoldChange <0)))))


  pet_v_gyn <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_PET_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2)
  pet_v_sta <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/STA_v_PET_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05 ) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  pet_v_sep <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/PET_v_SEP_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  pet_info <- data.frame('Comparison'=c('PET vs GYN','PET vs STA', 'PET v SEP'), 
  'Genes_Up'=c(nrow(pet_v_gyn %>% filter(log2FoldChange < 0)),nrow(pet_v_sta %>% filter(log2FoldChange < 0)),nrow(pet_v_sep %>% filter(log2FoldChange >0))),
  'Genes_Down'=c(-1*(nrow(sta_v_gyn %>% filter(log2FoldChange >0))),-1*(nrow(pet_v_sta %>% filter(log2FoldChange >0 ))),-1*(nrow(pet_v_sep %>% filter(log2FoldChange <0)))))

   sep_v_gyn <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/GYN_v_SEP_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2)
  sep_v_sta <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/STA_v_SEP_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05 ) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  sep_v_pet <- read.csv('~/Salba_RNA/results/organs/tissue/stage_8-9/PET_v_SEP_gtf_redo_no_lrt') %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  sep_info <- data.frame('Comparison'=c('SEP vs GYN','SEP vs STA', 'SEP v PET'), 
  'Genes_Up'=c(nrow(pet_v_gyn %>% filter(log2FoldChange < 0)),nrow(pet_v_sta %>% filter(log2FoldChange < 0)),nrow(pet_v_sep %>% filter(log2FoldChange <0))),
  'Genes_Down'=c(-1*(nrow(sta_v_gyn %>% filter(log2FoldChange >0))),-1*(nrow(pet_v_sta %>% filter(log2FoldChange > 0 ))),-1*(nrow(pet_v_sep %>% filter(log2FoldChange > 0)))))

 
  gyn_p <- ggplot(gyn_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    ylim(-5000, 5000)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Gynoecium')+
    ylab("Gene count")

  sta_p <- ggplot(sta_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    ylim(-5000, 5000)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Stamen')+
    ylab("Gene count")

    pet_p <- ggplot(pet_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    ylim(-5000, 5000)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Petal')+
    ylab("Gene count")  

  sep_p <- ggplot(gyn_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    ylim(-5000, 5000)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    )+ theme+
    ggtitle('Sepal')+
    ylab("Gene count")

  p <- ggpubr::ggarrange(gyn_p, sta_p, pet_p, sep_p, nrow = 1, ncol=4)
  ggsave('~/Salba_RNA/results/organs/tissue/stage_8-9/DEGs_lollipop.pdf', p,width = 210, height = 150, units = "mm", limitsize = F, dpi=500)


# BUD DEG lollipop ----
  six_v_eight <- read.csv('~/Salba_RNA/results/buds/stage_6_v_stage_8_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  six_v_ten<- read.csv('~/Salba_RNA/results/buds/stage_6_v_stage_10_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
  six_info <- data.frame('Comparison'=c('Stage 6-7 vs Stage 8-9','Stage 6-7 vs Stage 10-11'), 
  'Genes_Up'=c(nrow(six_v_eight %>% filter(log2FoldChange >0)),nrow(six_v_ten %>% filter(log2FoldChange >0))),
  'Genes_Down'=c(-1*(nrow(six_v_eight %>% filter(log2FoldChange <0))),-1*(nrow(six_v_ten %>% filter(log2FoldChange <0)))))


   eight_v_six <- read.csv('~/Salba_RNA/results/buds/stage_6_v_stage_8_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
   eight_v_ten <- read.csv('~/Salba_RNA/results/buds/stage_8_v_stage_10_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
   eight_info <- data.frame('Comparison'=c('Stage 8-9 vs Stage 6-7','Stage 8-9 vs Stage 10-11'), 
   'Genes_Up'=c(nrow(six_v_eight %>% filter(log2FoldChange < 0)),nrow(six_v_ten %>% filter(log2FoldChange >0))),
   'Genes_Down'=c(-1*(nrow(six_v_eight %>% filter(log2FoldChange >0))),-1*(nrow(six_v_ten %>% filter(log2FoldChange <0)))))


   ten_v_six <- read.csv('~/Salba_RNA/results/buds/stage_6_v_stage_10_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
   ten_v_eight <- read.csv('~/Salba_RNA/results/buds/stage_8_v_stage_10_gtf_redo_no_lrt') %>%
   dplyr::select(gene_id, log2FoldChange, padj) %>% na.omit %>% filter(padj < 0.05) %>% filter(log2FoldChange >2 |log2FoldChange < -2 )
   ten_info <- data.frame('Comparison'=c('Stage 10-11 vs Stage 6-7','Stage 10-11 vs Stage 8-9'), 
   'Genes_Up'=c(nrow(six_v_eight %>% filter(log2FoldChange < 0)),nrow(six_v_ten %>% filter(log2FoldChange >0))),
   'Genes_Down'=c(-1*(nrow(six_v_eight %>% filter(log2FoldChange >0))),-1*(nrow(six_v_ten %>% filter(log2FoldChange <0)))))


  six_p <- ggplot(six_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Bud Stage 6-7')+
    ylim(-3500, 3500)+
    ylab("Gene count")

  eight_p <- ggplot(eight_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Bud Stage 8-9')+
    ylim(-3500, 3500)+
    ylab("Gene count")

  ten_p <- ggplot(ten_info) +
    geom_segment( aes(x=Comparison, xend=Comparison, y=Genes_Up, yend=Genes_Down), color="grey") +
    geom_point( aes(x=Comparison, y=Genes_Up), color="#40B0A6",size=5 ) +
    geom_point( aes(x=Comparison, y=Genes_Down), color="#E1BE6A", size=5 ) +
    geom_hline(yintercept = 0,linetype = "dashed")+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle=90),
      axis.title.y = element_blank()
    ) + theme+
    ggtitle('Bud Stage 10-11')+
    ylim(-3500, 3500)+
    ylab("Gene count")  




ggsave('~/Desktop/test_plot.pdf', eight_p)




  p <- ggpubr::ggarrange(six_p, eight_p, ten_p, nrow = 1, ncol=3)
  ggsave('~/Salba_RNA/results/buds/DEGs_lollipop.pdf', p,width = 210, height = 150, units = "mm", limitsize = F, dpi=500)
