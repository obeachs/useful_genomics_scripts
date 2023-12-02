library(dplyr)
library(ggplot2)
library(tidyverse)

d <- read.csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/disorder_stats/AT_TRY_IUPRED.tsv', sep='\t')



palette("Paired")

# GGplot2 theme
theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  #plot.margin = unit(c(1, 1, 1, 1), "line")
)
prep_for_plot <- function(data){
    outname <- paste0(tools::file_path_sans_ext(data),".png")
    df <- read.csv(data,sep='\t')
    print(names(df))
    print(head(df))
    columns <- c("Measure","disScore","AA","POS") 
    temp_1 <- data.frame(matrix(nrow = length(df$X..POS), ncol = length(columns)))
    colnames(temp_1) <- columns
    temp_1 <- dplyr::mutate(temp_1, Measure='Anchor2', disScore=df$ANCHOR.SCORE, AA=df$AMINO.ACID, POS=df$X..POS)
    temp_2 = data.frame(matrix(nrow = length(df$X..POS), ncol = length(columns)))
    colnames(temp_2) = columns
    temp_2 <- dplyr::mutate(temp_2, Measure='IUPred', disScore=df$IUPRED.SCORE, AA=df$AMINO.ACID, POS=df$X..POS)
    new_df <- rbind(temp_1,temp_2)
    print(head(new_df))
    p <- ggplot(data=new_df, aes(x=POS, y=disScore, color=Measure) ) +
    geom_line(linewidth=3) +
    geom_hline(linetype = "dashed",yintercept=1.0)+
    geom_hline(yintercept=0.5)+
    geom_hline(linetype = "dashed",yintercept=0.0) + 
    xlab('Amino Acid Position')+
    ylab('Disorder Score')+
    scale_x_continuous(limits = c(0,max(new_df$POS)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,1.1), expand = c(0, 0)) +
    theme
    ggsave(outname,p, height=10, width=20, )
}


files <- list.files('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/disorder_stats/',full.names=T, pattern='.tsv')
lapply(files, prep_for_plot)
    

p <- ggplot(data=n, aes(x=POS, y=disScore, color=Measure)) +
    geom_line(linewidth=3) +
    geom_hline(linetype = "dashed",yintercept=1.0)+
    geom_hline(yintercept=0.5)+
    geom_hline(linetype = "dashed",yintercept=0.0) + 
    xlab('Amino Acid Position')+
    scale_x_continuous(limits = c(0,max(n$POS)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,1.1), expand = c(0, 0)) +
    theme

ggsave('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/disorder_stats/AT_TRY_IUPRED.png',p, height=10, width=20, )
