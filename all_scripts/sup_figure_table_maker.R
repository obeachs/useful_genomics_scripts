  source('/Volumes/sesame/joerecovery/scripts/formattable_functions_and_tables.R')
library(ggplot2)
library(cowplot)

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


nb.cols <- 8
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(8)
g <- tibble('Chromosome'=samps, 'Count'=count) %>% mutate(ratio=round((count/sum(count)*100)))
write.csv(g,'~/thesis_figs_and_tables/sup/weigel_hic_interactions.csv', row.names=F, quote=F)
make_nice_table(g,outname = '~/thesis_figs_and_tables/sup/weigel_hic_interactions.png')
raw <- ggplot(g, aes(x = 2, y = Count, fill = Chromosome)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = mycolors) +
  theme_void()+
  xlim(0.5, 2.5)+
  annotate(geom = 'text', x = 0.5, y = 0,size=10, label = paste0('Total: ',sum(g$Count)))+
  theme(legend.text = element_text(size=30))+
  theme(legend.title = element_text(size=30))+
  geom_text(aes(label = paste0(ratio,'%')),position = position_stack(vjust = 0.5),size=5, color = "white") +
  theme(axis.text = element_text(size = 5, colour="white"))


ggsave('~/thesis_figs_and_tables/sup/weigel_hic_interactions_donut.pdf', raw, height=10, width=10)

hicup <- ggplot(g, aes(x = 2, y = hicup_count, fill = Chromosome)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  scale_fill_manual(values = mycolors) +
  theme_void()+
  xlim(0.5, 2.5)+
  annotate(geom = 'text', x = 0.5, y = 0,size=10, label = paste0('Total: ',sum(g$HiCUP_Count)))+
   theme(legend.text = element_text(size=30))+
   theme(legend.title = element_text(size=30))


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




