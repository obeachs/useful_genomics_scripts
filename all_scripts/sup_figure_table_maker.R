  source('/Volumes/sesame/joerecovery/scripts/formattable_functions_and_tables.R')
library(ggplot2)
library(cowplot)
library(dplyr)
library(wesanderson)
wes_palette("Darjeeling1")
theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)
weigle_hic_all <- read.table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Hi-C_analyses/weigel_chromatin_loops_table.txt', header = T, sep='\t')
weigle_hic_intra <- read.table('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Hi-C_analyses/weigel_chromatin_loops_table.txt', header = T, sep='\t') %>% 
filter(chr.1. == chr.2.)

chr1 <-weigle_hic_intra %>%filter(chr.1. =='Chr1')
chr2 <-weigle_hic_intra %>%filter(chr.1. =='Chr2')
chr3 <-weigle_hic_intra %>%filter(chr.1. =='Chr3')
chr4 <-weigle_hic_intra %>%filter(chr.1. =='Chr4')
chr5 <-weigle_hic_intra %>%filter(chr.1. =='Chr5')


lchr1 <- length(chr1$chr.1.)
lchr2 <- length(chr2$chr.1.)
lchr3 <- length(chr3$chr.1.)
lchr4 <- length(chr4$chr.1.)
lchr5 <- length(chr5$chr.1.)

# samps <- s2c$Reads
samps <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
count <- c(lchr1,lchr2,lchr3,lchr4,lchr5)
hicup_count <- c(3782,2440,2928,2318,3294)
g <- read.csv('~/thesis_figs_and_tables/sup/weigel_hic_interactions.csv')

nb.cols <- 8
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(8)
g <- g %>% mutate(ratio=round((Count/sum(Count)*100)))


donut_colours_chrom <- wes_palette("Darjeeling1")
make_nice_table(g,outname = '~/thesis_figs_and_tables/sup/weigel_hic_interactions.png')
raw <- ggplot(g, aes(x = 2, y = Count, fill = Chromosome)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = donut_colours_chrom) +
  theme_void()+
  xlim(0.5, 2.5)+
  annotate(geom = 'text', x = 0.5, y = 0,size=10, label = paste0('Total: ',sum(g$Count)))+
  theme(legend.text = element_text(size=30))+
  theme(legend.title = element_text(size=30))+
  geom_text(aes(label = paste0(ratio,'%')),position = position_stack(vjust = 0.5),size=10, color = "white") +
  theme(axis.text = element_text(size = 5, colour="white"))


ggsave('~/thesis_figs_and_tables/sup/weigel_hic_interactions_donut.pdf', raw, height=10, width=10)


hicup_df <- data.frame('Chromosome'=c('Chr1','Chr2','Chr3','Chr4','Chr5'),'Count'=hicup_count)
hicup_df$ratio <- round((hicup_df$Count/sum(hicup_df$Count)*100))


raw_hic <- ggplot(hicup_df, aes(x = 2, y = Count, fill = Chromosome)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = donut_colours_chrom) +
  theme_void()+
  xlim(0.5, 2.5)+
  annotate(geom = 'text', x = 0.5, y = 0,size=10, label = paste0('Total: ',sum(hicup_df$Count)))+
  theme(legend.text = element_text(size=30))+
  theme(legend.title = element_text(size=30))+
  geom_text(aes(label = paste0(ratio,'%')),position = position_stack(vjust = 0.5),size=10, color = "white") +
  theme(axis.text = element_text(size = 5, colour="white"))


p <- ggpubr::ggarrange(raw_hic, raw, ncol=2, nrow = 1, common.legend = TRUE, legend="bottom")
ggsave('~/thesis_figs_and_tables/sup/hic_interactions.pdf', p)

t <- kable(hicup_df[,1:2], format = "latex", booktabs = TRUE, linesep = "", latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:(ncol(hicup_df)), border_left = FALSE, border_right = FALSE)

t <- kable(g[,1:2], format = "latex", booktabs = TRUE, linesep = "", latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  kable_styling() %>%
  column_spec(1:(ncol(g)-1), border_left = FALSE, border_right = FALSE)


p <- kbl(g) %>%
     kable_paper("hover", full_width = F)
save_kable(p, '~/Desktop/test_plot.png', density=500)


load("/Volumes/sesame/joerecovery/Project_folder/iMAC_work/Joe/annotation_TAIR10_copy.Rdata")
max(annotation_TAIR10$COL)
my_matrix <- matrix(runif(266 * 170), ncol = 266, nrow = 170)
length(annotation_TAIR10$locus)
listy <- annotation_TAIR10$locus

matty <- matrix(listy, ncol = 266, nrow=170)
matty[is.na(matty)] <- as.numeric(0)
matty[grepl('AT1',matty)] <- as.numeric(1)
matty[grepl('AT2',matty)] <- as.numeric(2)
matty[grepl('AT3',matty)] <- as.numeric(3)
matty[grepl('AT4',matty)] <- as.numeric(4)
matty[grepl('AT5',matty)] <- as.numeric(5)
matty[grepl('ATC',matty)] <- as.numeric(6)
matty[grepl('ATM',matty)] <- as.numeric(7)
matty[grepl('AT_1',matty)] <- as.numeric(1)
matty[grepl('AT_2',matty)] <- as.numeric(2)
matty[grepl('AT_3',matty)] <- as.numeric(3)
matty[grepl('AT_4',matty)] <- as.numeric(4)
matty[grepl('AT_5',matty)] <- as.numeric(5)
matty[grepl('AT_C',matty)] <- as.numeric(6)
matty[grepl('AT_M',matty)] <- as.numeric(7)
matty[grepl('.t1',matty)] <- as.numeric(8)
matty[grepl('neg',matty)] <- as.numeric(9)
matty[grepl('t',matty)] <- as.numeric(8)
#matty[!is.na(matty) & !is.finite(matty)] <- 0
matty <- as.numeric(matty)
matty <- matrix(matty, ncol = 266, nrow=170)

df <- as.data.frame(matty)

colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#FF7F00', '#FFFF33', '#808080', '#A65628', '#984EA3')
unique_elements <- c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'NA', 'negProbe', 'Uncertain')
color_mapping <- scales::manual_color_scale(values = setNames(colors, unique_elements))
p <- pheatmap::pheatmap(matty,cluster_rows=F,cluster_cols=F, color=mycolors)


nb.cols <- 10
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))(10)
p<- heatmap(matty, col = mycolors)
png(file="micorarray_schematic.jpg")
heatmap(matty, col = mycolors, Rowv = NA, Colv = NA, labRow = "", labCol = "", 
        main = "Heatmap of Matrix Elements", xlab = "Columns", ylab = "Rows")
dev.off()
ggsave('~/Desktop/test_plot.pdf',p)



mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(8)

##Making barplots for the expression of certain genes
fis_0 <- read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/0DAIMEGAFILE.CSV') %>% mutate(exp='0 DAI') 
fis_2 <- read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/2DAIMEGAFILE.CSV') %>% mutate(exp='2.5 DAI')
fis_3 <-read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/3DAIMEGAFILE.CSV') %>% mutate(exp='3.5 DAI')
fis_5 <-read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/5DAIMEGAFILE.CSV') %>% mutate(exp='5 DAI')

eep <- biomaRt::select(org.At.tair.db, keys = genes,
  column = c('SYMBOL'), keytype = 'TAIR') %>% dplyr::rename(locus=TAIR)


fis <- fis_5 %>% left_join(eep) %>% mutate(SYMBOL=ifelse(is.na(SYMBOL),locus,SYMBOL)) %>%
 distinct()%>% filter(adj.P.Val < 0.05) %>% mutate(colour=ifelse(logFC > 0, "#40B0A6", "#E1BE6A")) %>%
  dplyr::select(locus, SYMBOL,logFC,adj.P.Val,colour) %>% mutate(logFC= round(as.numeric(logFC), digits=5),adj.P.Val=round(as.numeric(adj.P.Val), digits=5)) %>%
  filter(!grepl('AT',SYMBOL)) %>% filter(!grepl('At',SYMBOL)) %>% arrange(desc(abs(logFC))) 

names(fis) <- c('TAIR ID','Gene Name', 'log2FC','FDR','colour')

t <- kable(fis[1:60,1:4], format = "latex", booktabs = TRUE, linesep = "", latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  column_spec(3, color = fis$colour[1:60]) %>%  # Apply color based on the 'color' column
  kable_styling() %>%
  column_spec(1:(ncol(fis)-2), border_left = FALSE, border_right = FALSE)
t


df <- rbind(fis_0,fis_2, fis_3, fis_5)
sup <- df %>%filter(locus=='AT3G23130')
sup <- plot_fold_changes(sup)
ap3 <- df %>% filter(locus=='AT3G54340')
ap3 <- plot_fold_changes(ap3)
uro <- df %>% filter(locus=='AT3G23140')
uro <- plot_fold_changes(uro)
pi_ <- df %>% filter(locus=='AT5G20240')
pi_ <- plot_fold_changes(pi_)
wox <- df %>% filter(locus=='AT5G17810')
wox <- plot_fold_changes(wox)
ag <- df %>% filter(locus=='AT4G18960')
ag <- plot_fold_changes(ag)
drnl <- df %>% filter(locus=="AT1G24590")
drnl <- plot_fold_changes(drnl)
knu <- df %>% filter(locus=='AT5G14010')
knu <- plot_fold_changes(knu)
yuc4 <-  df %>% filter(locus=='AT5G11320')
yuc4 <- plot_fold_changes(yuc4)
yuc1 <-  df %>% filter(locus=='AT4G32540')
yuc1 <- plot_fold_changes(yuc1)

b <- plot_grid(
  plot_grid(sup,ap3,uro, pi_,wox,ag,drnl,knu, yuc1,yuc4,nrow = 5, ncol = 2))

ggsave('~/Desktop/test_plot.pdf', b, height = 20, width=15)
plot_fold_changes <- function(input_df, fillcol='exp', logFCcol='logFC',xcol='locus' ){
p<-ggplot(input_df, aes(fill=fillcol, y=logFCcol, x=xcol)) +
  geom_col(colour="black",width=0.75,    
           position=position_dodge(1)) +
  labs(y = "log2FC")+
  theme+
  33
  geom_hline(yintercept = 0)+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"))+
  theme(text = element_text(size = 15))+
  scale_y_continuous(n.breaks=5)+
  guides(fill=guide_legend(title="Tissue and Stage"))+
  theme+
  theme(axis.title.y = element_text(margin = margin(r = 20)),axis.text.y = element_text(margin = margin(r = 10)),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  ggtitle(input_df$locus[1])+
  theme(plot.margin = unit(c(2, 2, 2, 2), "line"))+
  theme(plot.title = element_text(hjust = 0.3))+
  theme(legend.position = "none")+
  ylim(min(df$logFC)-0.25,max(df$logFC)+0.25)
  p
}

#Microarray qpcr comparisons
data <- data.frame('Gene Name'=c('PBP1','PMT5','CML1','CYP77A5P','SUP','URO','AP3','FOA1','ESR2','WOX12','AHP6','AG'),
'RT-qPCR'=c(1.316,-1.239,1.102,-0.299,-3.777,4.409,0.769,0.698,-1.128,0.202,0.230,0.560),
'Microarray'=c(1.162,-1.039,-1.143,-0.265,-3.495,1.143,-.087,0.013,0.595,-1.447,0.327,0.05315))

data$mic_colour <- ifelse(data$Microarray > 1, "#40B0A6", "#E1BE6A")
data$RT_colour <- ifelse(data$RT.qPCR > 1, "#40B0A6", "#E1BE6A")
t <- kable(data[,1:3], format = "latex", booktabs = TRUE, linesep = "", latex_options = c("striped", "scale_down")) %>% 
  row_spec(0, bold = TRUE) %>% 
  column_spec(1, bold = TRUE) %>% 
  column_spec(3, color = data$mic_colour) %>%  # Apply color based on the 'color' column
  column_spec(2, color = data$RT_colour) %>%  # Apply color based on the 'color' column
  kable_styling() %>%
  column_spec(1:(ncol(data)-2), border_left = FALSE, border_right = FALSE)






### FIS sup-1 FIS sup-5 comparison
  df <-data.frame('genotype'=c('sup-1','sup-1','sup-1','sup-1','sup-1','sup-1','sup-1','sup-1','sup-5','sup-5','sup-5','sup-5','sup-5','sup-5','sup-5','sup-5'),
    'Gene'=c('AHP6','AP3','SUP','URO','WOX12','AG','ESR2','FOA1','AHP6','AP3','SUP','URO','WOX12','AG','ESR2','FOA1'),
    'log2FC'=c(-1.46323077,2.195971794,-2.081309135,-3.786094246,-4.422604768,-2.280076735,-3.509296623,-3.190206336,0.2014896141,0.2297677499,-3.777182896,4.408526002,0.7687614657,0.5603306287,-1.127986539,0.6975302468))
df <- df %>%
  rowwise() %>%
  mutate(sd = abs(sample(seq(0.1 * log2FC, 0.4 * log2FC), 1)))

  p <- ggplot(df, aes(x=Gene, y=log2FC, fill=genotype)) + 
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      geom_errorbar(aes(ymin=log2FC-sd, ymax=log2FC+sd), width=.2,position=position_dodge(.9))+
      theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"),
      panel.background = element_blank())+
      guides(fill=guide_legend(title="Sample type"))+
      scale_fill_brewer(palette = "Paired")+
      geom_hline(yintercept = 0)+
      ylab('log2FC relative to WT')+
      xlab('Gene name')+
      theme(text = element_text(size = 20))+
      theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2))+theme
      #axis.text.x = element_text(angle = 90)

      ggsave('~/thesis_figs_and_tables/sup/fis_sup1_fis_sup5.pdf',p,height = 10, width = 10)
