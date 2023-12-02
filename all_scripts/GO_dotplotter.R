library(dplyr)
library(ggplot2)

#Loading data and selecting the columns we want and giving them sane column names
#Avoid messing with these next 6 or so lines
data <- read.csv('/Volumes/sesame/movers/DEG Enrichment data.csv')
data <- data[,c(1,6,7)]
names(data) <- c('GO_term','Enrichment','Raw_pvalue')



#From kegg
data <- as.data.frame(ego_rescues@result)%>% dplyr::select(Description, Count, pvalue)
names(data) <- c('GO_term','Enrichment','Raw_pvalue')








data$Enrichment <- as.numeric(data$Enrichment)
data <- dplyr::arrange(data, desc(Enrichment)) %>% dplyr::slice(1:50)
data$GO_term <- fct_rev(factor(data$GO_term, levels = data$GO_term))
data$GO_term <- 

#Change around the colours in the gradient if needed
p <- ggplot(data = data, aes(x = Enrichment, y = GO_term, 
                        color = Raw_pvalue, size = Enrichment)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("GO Terms - Biological process") + 
  xlab("Enrichment") + 
  ggtitle("GO enrichment analysis")


#Change name as you please
ggsave('~/coolplot.pdf', p)

