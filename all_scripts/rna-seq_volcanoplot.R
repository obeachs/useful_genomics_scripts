library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for the labelling
library(ggplot2)

#Simple script to get volcano plots from RNA-seq differential expression analysis data
#Keep an eye on the comments to see if there are any changes by the user needed
2

#Read in the data - best way to get filepath on mac is to go to finder and right click then hit 'Get Info' and copy and paste the shown path
data <- read.csv('/Users/josephbeegan//AC 7.7.23 BTZ 6HRS Raw Data.xlsx - data3_2_vs_1_fc2_&_raw.p.csv')
/Users/josephbeegan/Downloads/
#VERY IMPORTANT -Change the words in the quotation marks here to match the names of the columns in your dataset
#Make sure it's in the order of the equivalent of 1.GeneName 2.Fold Change 3. Pvalue
#The order of the columns in the dataset does not matter at all, just the order that they're written down in the function below
data <- select(data,'Gene_Symbol','X2.1.fc','X2.1.raw.pval') 
#Adding a new column to the data to say if the gene is upregulated or downreulgated


data <- data %>%
  mutate(Direction = case_when(
    data[,2]  >1 ~ "UP",
    data[,2] <1  ~ "DOWN"
    ))

#Adding a new column to the data for labelling of speciifc genes on the plot
#As an example just taking the top thirty here but if needed can do any GOI
data$labelme <- ifelse(data[,1] %in% head(data[order(data[,3]), "Gene_Symbol"], 30), data[,1], NA)
data$loggy <- abs(data[,2])
data$loggy <- log2(data$loggy)
data$loggy[data$Direction == 'DOWN'] <- -data$loggy[data$Direction == 'DOWN']
data$X2.1.raw.pval[data$X2.1.raw.pval==0] <- 0.0000000001

#This looks complicated but just plotting the data and adding some lines, if you wanted to mess around with the
#Look of the label, just change the numbers in the 'geom_text_repel'
#If you want tochange the tile of the plot, change the words in the quotations marks near ggtitle

outplot <- ggplot(data = data, aes(x = loggy, y = -log10(data[,3]), col = Direction, label=labelme)) +
     geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(position = 'jitter') +
   scale_color_manual(values = c("red", "#28bf21"), # to set the colours of our variable  
                     labels = c("Downregulated", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
 coord_cartesian(ylim = c(0, 12), xlim = c(-10, 15)) + # since some genes can have minuslog10padj of inf, we set these limits

  theme(plot.margin = margin(30, 30, 30, 30, "pt"))+
  labs(color = 'DEG status', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
       ylim(0,12) +
theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )+ 
  xlab("log2FC") +
  ylab("-log10pvalue") +
  ggtitle('UP vs DOWN')+ geom_label_repel(fill = "white",box.padding = 1,max.overlaps = Inf,arrow = arrow(length = unit(0.005, "npc"))) # To show all labels 




#This is saving the plot as a pdf - just change the path and name to whatever you want it to be
ggsave('~/Desktop/aoifeplot.pdf', outplot, height = 30, width = 10, limitsize = F)

