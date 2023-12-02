library(ggplot2)
library(dplyr)

theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)





df2 <- read.csv('~/Documents/Salba_qPCRs/JB_Salba_ORGANS_TTG1-305-306_qPCR_analysis_copy.txt',sep='\t')
df2$sam_stage <- sub("MID", "_stage_8-9", df2$sam_type)
df2$sam_stage <- sub("LATE", "_stage_10-11", df2$sam_stage)
df2$forx <- paste(df2$sam_name, df2$sam_stage,sep = '_')
df2$ordercol <- df2$forx
df2$ordercol <- sub("stage_8-9", "EARLY", df2$ordercol)
df2$ordercol <- sub("stage_10-11", "LATE", df2$ordercol)
df2$ordercol <- paste(df2$sam_name,df2$ordercol, sep = '_')
df2 <- arrange(df2, ordercol)
custom_order <- df2[order(df2$ordercol), "forx"]
df2$forx <- factor(df2$forx, levels = custom_order)
df2 <- filter(df2, cp < 1000)
deep <- df2
deep <- df2 %>% filter(grepl('GL1a_', forx))
deep <- df2 %>% filter(grepl('GL1a-2', forx))
Sal06g09840_GIS3 <- df2 %>% filter(grepl('GIS3_', forx))
deep <- df2 %>% filter(grepl('GL2', forx))
deep <- df2 %>% filter(grepl('Sal04g07070L', forx))
deep <- df2 %>% filter(grepl('GL1b-2', forx))
deep <- df2 %>% filter(grepl('Sal07g02330L-GIS-2',forx))
deep <- df2 %>% filter(grepl('Sal07g02330L-GIS_',forx))
deep <- df2 %>% filter(grepl('WER_',forx))
deep <- df2 %>% filter(grepl('WER-2',forx))
deep <- df2 %>% filter(grepl('Sal10g28990L-GL1b_',forx))
deep <- df2 %>% filter(grepl('Sal06g09840L-GIS3_',forx))
deep <- df2 %>% filter(grepl('Sal04g13960L-GL1a-2',forx))
deep <- df2 %>% filter(grepl('Sal04g13960L-GL1a_',forx))

Sal04g07070L_CML42 <- df2 %>% filter(grepl('Sal04g07070L', forx))
Sal08g05440L_XI1 <- df2 %>% filter(grepl('Sal08g05440L-XI1', forx))
scale_fill_brewer(palette = "Paired")
#df2$sam_type <- factor(df2$sam_type, levels = c("stage_6-7", "stage_8-9", "stage_10-11"))
#df2$forx <- paste(df2$sam, df2$sam_type,sep = '_')

p <- ggplot(Sal08g05440L_XI1, aes(x=forx, y=cp, fill=sam_stage)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=cp-sd, ymax=cp+sd), width=.2,position=position_dodge(.9))+
    theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"),
    panel.background = element_blank())+
    guides(fill=guide_legend(title="Sample type"))+
    scale_fill_brewer(palette = "Paired")+
    ylab('Value relative to reference TUB1')+
    xlab('Gene name and Sample')+
    theme(text = element_text(size = 20))+
    theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2))
    #axis.text.x = element_text(angle = 90)
print(p)
    ggsave('~/thesis_figs_and_tables/salba/Sal08g05440L_XI1_qpcr.pdf',p,height = 10, width = 10)

title <- substr(deep$ordercol[1],1,12)
deep$sam_stage <- gsub('stage','Stage',deep$sam_stage)
p <- ggplot(deep, aes(x=sam_stage, y=cp, fill=sam_stage)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=cp-sd, ymax=cp+sd), width=.2,position=position_dodge(.9))+
    theme(axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"),
    panel.background = element_blank())+
    guides(fill=guide_legend(title="Sample type"))+
    scale_fill_brewer(palette = "Paired")+
    ylab('Value relative to reference TUB1')+
    xlab('Sample')+
    theme+
    theme(text = element_text(size = 20))+
   scale_x_discrete(limits = c(
    "GYN_Stage_8-9", "GYN_Stage_10-11",
    "PET_Stage_8-9", "PET_Stage_10-11",
    "SEP_Stage_8-9", "SEP_Stage_10-11",
    "STA_Stage_8-9", "STA_Stage_10-11"))+
    guides(fill=guide_legend(title="Tissue and Stage"))+
  ggtitle(paste("RT-qPCR:", title, sep=' '))
ggsave('~/thesis_figs_and_tables/salba/Sal04g18500-TTG1_qpcr.pdf',p,height = 10, width = 10)







df <- read.csv('~/Documents/Salba_qPCRs/JB Salba 27-4-23 2_qPCR_analysis.tsv', sep = '\t') 
q <- ggplot(df2, aes(x=sam, y=cp, fill=sam_type)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=cp-sd, ymax=cp+sd), width=.2,position=position_dodge(.9))+
    theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90))+
    scale_fill_brewer(palette = "Paired")+
    ylab('Value relative to reference TUB1')+
    xlab('Sample name')
    #axis.text.x = element_text(angle = 90)
print(p)ls -l
ggpubr::ggarrange(p,q, ncol = 1, nrow = 2)
# Finished bar plot
p+labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))
df2



#Chip_qPCR plots

IP2	EXP 	2.345669898
IP3	EXP 	1.729074463
IP4	EXP 	1.140763716
"URO 394/395",
"URO 396/397",
"YUC1 384/385",
"YUC4 392/393",
"AP3 514/515",
"1.69"
"1.07"
"0.208"
"0.2624"
"0.238"
"1.31",
"1.24",
"0.95",
"0.76",
"0.67",
"1.4",
"1.18",
"1.32",
"1.04",
"0.9"





df <- data.frame (IP  = c("IP2", "IP3", "IP4","IP2", "IP3", "IP4","IP2", "IP3", "IP4","IP2", "IP3", "IP4","IP2", "IP3", "IP4"),
                  Gene = c("URO 394/395",
"URO 396/397",
"YUC1 384/385",
"YUC4 392/393",
"AP3 514/515",
"URO 394/395",
"URO 396/397",
"YUC1 384/385",
"YUC4 392/393",
"AP3 514/515",
"URO 394/395",
"URO 396/397",
"YUC1 384/385",
"YUC4 392/393",
"AP3 514/515"),
                  foldchange =c("1.69",
"1.07",
"0.208",
"0.2624",
"0.238",
"1.31",
"1.24",
"0.95",
"0.76",
"0.67",
"1.4",
"1.18",
"1.32",
"1.04",
"0.9")
                  )
q <- ggplot(df, aes(x=Gene, y=as.numeric(foldchange), fill=IP)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90))+
    scale_fill_brewer(palette = "Paired")+ 
    theme+
    ylab('Value relative to Input')+
    xlab('Gene name') +
    geom_hline(yintercept = 1)
