library(ggplot2)
library(dplyr)


df2 <- read.csv('~/Documents/Salba_qPCRs/JB Salba TTG1 organs_qPCR_analysis.tsv', sep = '\t') 
df2$forx <- paste(df2$sam, df2$sam_type,sep = '_')


p <- ggplot(df2, aes(x=forx, y=cp, fill=sam_type)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=cp-sd, ymax=cp+sd), width=.2,position=position_dodge(.9))+
    theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"),
    panel.background = element_blank())+
    guides(fill=guide_legend(title="Sample type"))+
    scale_fill_brewer(palette = "Paired")+
    ylab('Value relative to reference TUB1')+
    xlab('Gene name and Sample')
    #axis.text.x = element_text(angle = 90)
print(p)
ggsave('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/final_figures/qpcr/ttg1_mid_organs.pdf',p)

df <- read.csv('~/Documents/Salba_qPCRs/JB Salba 27-4-23 2_qPCR_analysis.tsv', sep = '\t') 
q <- ggplot(df, aes(x=sam, y=cp, fill=sam_type)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=cp-sd, ymax=cp+sd), width=.2,position=position_dodge(.9))+
    theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90))+
    scale_fill_brewer(palette = "Paired")+
    ylab('Value relative to reference TUB1')+
    xlab('Sample name')
    #axis.text.x = element_text(angle = 90)
print(p)
ggpubr::ggarrange(p,q, ncol = 1, nrow = 2)
# Finished bar plot
p+labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))
df2
