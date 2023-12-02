library(ggplot2)
library(dplyr)
library(tidyverse)
library(gggenes)
library(RColorBrewer)
library(pals)
library(drawProteins)
library(ggpubr)

library(wesanderson)
files <- list.files('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/', pattern = 'hitdata.txt', full.names = TRUE)
list_of_domains <- list()
for (file in files){
  data <- read.csv(file,skip = 7, sep = '\t', header = TRUE)
  data <- data[(data$Hit.type == 'specific' | data$Hit.type == 'superfamily'),]
  data <- dplyr::select(data, Query,Short.name, From, To)
  data <- dplyr::rename(data, Name = Short.name, Start=From, End=To)
  data$Query <- stringr::str_sub(data$Query,8,-1)
  list <- append(list,data$Name)
}
list <- unique(list)
colors <- rainbow(length(list))
custom <- data.frame(names = list, 
                         color = colors, 
                         stringsAsFactors = FALSE)
custom$names <- sample(custom$names)
custom
my_pal <- colorRampPalette(rev(brewer.pal(n = (length(list)), name = "PuRd")))
my_cutstom <- data.frame(names = list, 
                         color = my_pal(length(list)), 
                         stringsAsFactors = FALSE)

my_cutstom[nrow(my_cutstom) + 1,] = c("space","white")
pals::pal.bands(my_cutstom$color)

print_protein_domains <- function(hitdata,prot_len){
  name <- hitdata
  data <- read.csv(hitdata,skip = 7, sep = '\t', header = TRUE)
  data <- data[(data$Hit.type == 'specific' | data$Hit.type == 'superfamily'),]
  data <- dplyr::select(data, Query,Short.name, From, To)
  data <- dplyr::rename(data, Name = Short.name, Start=From, End=To)
  data$Query <- stringr::str_sub(data$Query,8,-1)
  data$width <-(data$End - data$Start)
  data <- arrange(data,width)
  print(data)
  outname <- str_sub(name,1,-5)
  outname <- paste(outname,'_figure.pdf', sep = '')
  p <-ggplot(data = data, aes(x = Query, y=prot_len)) + geom_bar(stat="identity", fill="grey", width = +0.05) +coord_flip()+labs(x="Gene Name", y="Position")
  #p <- p + geom_rect(data=data, mapping=aes(xmin=Start, xmax=End, ymin=Query, ymax=Query, fill=Name), alpha=0.5)
  p <-p +geom_segment(aes(x=Query, xend=Query, y=Start,z=width, yend=End, color=Name, linewidth=10,linetype='dashed', size=Name),position='stack')+ scale_colour_manual(name=my_cutstom$names, values = custom$color)+
  scale_fill_manual(scale_colour_manual(name=my_cutstom$names, values = custom$color,1))
  #scale_color_identity()
  # <-p +geom_segment(data=data, aes(x=Query, xend=Query, y=Start, yend=End, color=Name),position='stack', alpha = 0.5) # nolint: line_length_linter.
  #p <- annotate(geom = "rect", xmin = 1, xmax = 30,ymin = 0, ymax = 5,color = Name,fill = Name, alpha = 0.5)
  ggsave(outname, p)
}
lapply(files,print_protein_domains,100)

print_protein_domains_gggenes <- function(hitdata){
  name <- hitdata
  data <- read.csv(hitdata,skip = 7, sep = '\t', header = TRUE)
  data <- data[(data$Hit.type == 'specific' | data$Hit.type == 'superfamily'),]
  print(data)
  data <- dplyr::select(data, Query,Short.name, From, To)
  data <- dplyr::rename(data, Name = Short.name, Start=From, End=To)
  data$Query <- stringr::str_sub(data$Query,8,-1)
  data$width <- data$End - data$Start
  data <- arrange(data,desc(width))
  outname <- str_sub(name,1,-5)
  outname <- paste(outname,'_figure_gggenes.pdf', sep = '')
  p <- ggplot(data, aes(xmin = Start, xmax = End, y = Query, fill = Name)) +
  geom_gene_arrow() +
  #facet_wrap(~ Query, scales = "fixed", ncol = 1) +
  scale_fill_brewer(palette = "Set3") + theme_genes()
  ggsave(outname, p)
}
lapply(files,print_protein_domains_gggenes, 500)


beep <- read.csv(files[1], skip=7, sep = '\t', header = TRUE)
data <- dplyr::select(beep, Query,Short.name, From, To)
data <- dplyr::rename(data, Name = Short.name, Start=From, End=To)
data$Query <- stringr::str_sub(data$Query,8,-1)
data$width <- data$End - data$Start
data <- arrange(data,desc(width))
p <-ggplot(data = data, aes(x = Start, y=Query)) + geom_bar(stat="identity", fill="grey", width = +0.05) + geom_rect(data = data, aes(xmin = Start, ymin = End, xmax = Start + 10, ymax = End + 0.2,fill = Name), position = 'stack') # nolint
ggplot(data, aes(xmin = Start, xmax = End, y = Query, fill = Name)) +
  geom_gene_arrow() +
  #facet_wrap(~ Query, scales = "fixed", ncol = 1) +
  scale_fill_brewer(palette = "Set3") + theme_genes()



features <- list(
  list(name = "SANT", start = 73, end = 116, color = "red"),
  list(name = "Myb_DNA-binding", start = 71, end = 116, color = "blue"),
  list(name = "SANT superfamily", start = 71, end = 116, color = "green"),
  list(name = "REB1 superfamily", start = 17, end = 113, color = "orange"),
  list(name = "PLN03091 superfamily", start = 17, end = 125, color = "purple"),
  list(name = "PLN03212 superfamily", start = 17, end = 137, color = "pink")
)
drawProteins::draw_canvas(features) + 
  draw_protein_legend(position = "bottom")
rel <- drawProteins::get_features("Q8GV05")
drawProteins::get_features("Q8GV05") ->
    rel_json
drawProteins::feature_to_dataframe(rel_json) -> rel_data


pro <- read.csv('/Volumes/sesame/movers/iprscan5-R20230503-090833-0292-5219440-p2m.tsv', sep = '\t', header = F)
pro <- pro[c(1,3,7,8,13)]
names(pro) <- c('Gene','length','start','end','name')
pro <- pro[(pro$name != '-'),]
pro_alba <- pro[(pro$Gene == 'Sal02g14170L'),]
pro_alba$name_again <- pro_alba$name
pro_alba <- pro_alba[order(pro_alba$start, decreasing = F),]


# pro_alba$tier <- 1
# for(i in 1:(nrow(pro_alba)-1)) {
#   for(j in (i+1):nrow(pro_alba)) {
#     if(pro_alba$name[i] != pro_alba$name[j]) {
#       if(pro_alba$start[i] <= pro_alba$end[j] & pro_alba$start[j] <= pro_alba$end[i]) {
#         pro_alba$tier[i] <- pro_alba$tier[]+1
#       }
#     }
#   }
# }
# Convert the name column to numeric


pro_alba$start <- as.numeric(pro_alba$start)
pro_alba$end <- as.numeric(pro_alba$end)
# Sort the data frame in descending order based on the length of each domain
pro_alba <- pro_alba[order(pro_alba$start, decreasing = F),]
pro_alba <- pro_alba[order(pro_alba$end - pro_alba$start, decreasing = TRUE),]
pro_alba$name <- as.numeric(as.factor(pro_alba$name))
pro_alba$name <- pro_alba$name + 1
pro_alba[nrow(pro_alba) + 1,] = c(pro_alba$Gene[1],as.numeric(pro_alba$length[1]),as.numeric(0),as.numeric(pro_alba$length[1]),as.numeric(1),'Protein')
pro_alba$start <- as.numeric(pro_alba$start)
pro_alba$end <- as.numeric(pro_alba$end)
pro_alba$name <- as.numeric(pro_alba$name)
# Create a horizontal stacked bar chart
myguide <- guide_legend(keywidth = unit(3, "cm"))
ggplot(pro_alba, aes(y = name, x = end - start, fill = as.factor(name_again))) +coord_flip()+
  geom_rect(aes(xmin = start, xmax = end, ymin = name - 0.5, ymax =name + 0.5), 
            color = "black", size = 0.5) +
  scale_y_continuous(breaks = unique(pro_alba$name), labels = levels(pro_alba$name)) +
  scale_fill_discrete() +
  theme_bw() +
  theme(legend.text = element_text(size=10),legend.key.height = unit(5, 'cm'))+
  xlab("Domain length") +
  guides(fill=guide_legend(title="Domain name"))+
  ylab("") +
  ggtitle(paste("Domain structure of ", pro_alba$Gene[1]))



interpro <- list.files('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/interpro', full.names = T)

print(file)
pro <- read.csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/yang_assemblies/sinapis_alba/trichome/prot/interpro/Sal02g14170L_interpro.tsv', sep = '\t', header = F)
pro <- pro[c(1,3,7,8,13)]
names(pro) <- c('Gene','length','start','end','name')

pro <- pro[(pro$name != '-'),]
plot_list <- vector('list', length = interpro)
for (i in 1:length(unique(pro$Gene))){
  pro_alba <- pro[(pro$Gene == unique(pro$Gene)[i]),]
  print(pro_alba)
  pro_alba$name_again <- pro_alba$name
  pro_alba <- pro_alba[order(pro_alba$start, decreasing = F),]
  pro_alba$name <- as.numeric(as.factor(pro_alba$name))
  pro_alba$name <- pro_alba$name + 1
  pro_alba[nrow(pro_alba) + 1,] = c(pro_alba$Gene[1],as.numeric(pro_alba$length[1]),as.numeric(0),as.numeric(pro_alba$length[1]),as.numeric(1),'Protein') # nolint
  pro_alba$start <- as.numeric(pro_alba$start)
  pro_alba$end <- as.numeric(pro_alba$end)
  pro_alba$name <- as.numeric(pro_alba$name)
  # Create a horizontal stacked bar chart
  myguide <- guide_legend(keywidth = unit(3, "cm"))
  p <- ggplot(pro_alba, aes(y = name, x = end - start, fill = as.factor(name_again))) +coord_flip()+
    geom_rect(aes(xmin = start, xmax = end, ymin = name - 0.5, ymax =name + 0.5), 
              color = "black", size = 0.5) +
    scale_y_continuous(breaks = unique(pro_alba$name), labels = levels(pro_alba$name)) +
    scale_fill_discrete() +
    theme_bw() +
    theme(legend.text = element_text(size=10),legend.key.height = unit(2.5, 'cm'))+
    xlab("Protein length") +
    guides(fill=guide_legend(title="Domain name"))+
    ylab("") +
    scale_x_continuous(n.breaks=25) +
    scale_fill_brewer(palette = "Paired")+
    rotate() +
    scale_y_reverse() +
    theme(axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) + 
    ggtitle(paste("Domain structure of ", pro_alba$Gene[1]))
    out_name <- paste('~/Desktop/protein_out_', pro_alba$Gene[1], sep = '')
    out_name <- paste(out_name,'.pdf', sep = '')
    ggsave(out_name,p)
    plot_list[[i]] <- p
  }

ggsave('~/Desktop/protein_out.pdf',ggpubr::ggarrange(plotlist=plot_list, widths = c(3,3)))
