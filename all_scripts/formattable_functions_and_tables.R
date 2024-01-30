library(HELP)
library(dplyr)
library(Ringo)
library(org.At.tair.db)
library(biomaRt)
library(dplyr)
library(htmltools)
library(formattable)
library(webshot)
library(tidyr)
library(kableExtra)
#YOU NEED GOOGLE CHROME INSTALLED FOR THIS
#Function to export the table as png
export_formattable <- function(f, file, width = "50%", height = NULL, 
                               background = "white", delay = 0.2, vwidth=2000, vheight=750){
      w <- as.htmlwidget(f, width = width, height = height)
      path <- html_print(w, background = background, viewer = NULL)
      url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
      webshot2::webshot(url,
              file = file,
              selector = ".formattable_widget",
              delay = delay,
              vwidth = vwidth, vheight = vheight)
    }

italicize <- function(x) {
  formattable::color_tile("transparent", "transparent", 
                          x, format = "html")
}
#Function to make any values above 0 green and below 0 red
#Usage like this:
#formattable::formattable(df, list(
#  targetColName = color_formatter,
#   = formattable::color_tile("transparent","lightgreen")))
#Can be used on multiple columns
color_tile_mean <- function (...) {
  formatter("span", style = function(x) {
    style(display = "block",
          padding = "0 4px", 
          `border-radius` = "4px", 
          `background-color` = ifelse(x < 0 , '#f96955', '#aff581')) # Remember to change the colors!
  })}

#These two lines will make the titles of each column bold

format_columns_by_indices <- function(df, indices=c()) {
  if(length(indices) >0){
  formatted_df <- df
  for (col_index in indices) {
    formatted_df[, col_index] <- color_tile_mean()(formatted_df[, col_index])
  }
  formattable(formatted_df, list(indices),align = c("c", "c", "c"))
  }
  else{
  p <- formattable(df, align=c("c", "c", "c"),list(
    col.names = as.character(indices)))
  }
}
setwd('~/Salba_RNA/')




make_nice_table <- function(df,colnums=c(),outname,vwidth=2000, vheight=750){
  p <- format_columns_by_indices(df, indices = colnums)
  export_formattable(p, outname, vwidth, vheight)
}
# f <- formattable(mtcars)
# export_formattable(f,'~/Desktop/test_plot.pdf')

# formattable(mtcars) %>% 
#   as.htmlwidget() %>% 
#   htmlwidgets::saveWidget(file="test.html")
# #make_nice_table(df=df, colnums = c(2,3), outname ='/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/final_figures_and_tables/RT-qPCR_microarray_value_comparison.png')

  # alp_data <- read.csv('/Volumes/sesame/ALP_Omics/sample_ids.csv') %>% dplyr::select(Genotype, Condition, ExperimentType, Library) %>% distinct()
# alp_data <- alp_data[!grepl('Input', alp_data$Condition),]
# alp_data$Condition <- ifelse(grepl('H3', alp_data$Condition),alp_data$Condition, '') 
# alp_data$ExperimentType <- paste(alp_data$Condition, alp_data$ExperimentType,sep = ' ')
# alp_data <- alp_data %>% dplyr::select(Genotype, ExperimentType, Library) %>% dplyr::rename(Experiment=ExperimentType, 'Library type'= Library)
# row.names(alp_data) <- NULL



# # custom_formatter <- function(x) {
# #   formattable::as_paragraph(
# #     ifelse(x != "Col-0", formattable::as_css("font-style: italic;"), ""),
# #     x
# #   )
# # }
# alp_data$Genotype <- formattable::formattable(alp_data$Genotype, formatter = custom_formatter)
# formattable_table <- formattable::formattable(alp_data)


# make_nice_table(alp_data,outname = '~/Desktop/test_plot.pdf')
# formattable(alp_data)

# #Organ-specific floral genes
make_nice_table(data.frame(Organ=c('Gynoecium','Stamen','Petal','Sepal'),Count=c(1173,4099,469,1323)), outname = './results/flower/organ_specific_genes/stage_10-11_table.pdf')
make_nice_table(data.frame(Organ=c('Gynoecium','Stamen','Petal','Sepal'),Count=c(522,3273,504,2742)), outname = './results/flower/organ_specific_genes/stage_8-9_table.pdf')
make_nice_table(data.frame(Organ=c('Gynoecium','Stamen','Petal','Sepal'),Count=c(1251,5359,737,3088)), outname = './results/flower/organ_specific_genes/all_stages_table.pdf')

stage_8 <- data.frame(Organ=c('Gynoecium','Stamen','Petal','Sepal'),Count=c(1173,4099,469,1323))
stage_10 <- data.frame(Organ=c('Gynoecium','Stamen','Petal','Sepal'),Count=c(522,3273,504,2742))






# # df <- read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/final_figures_and_tables/RT-qPCR_microarray_value_comparison.csv')


# # format_columns_by_indices <- function(df, indices=c()) {
# #   if(length(indices) >0){
# #   formatted_df <- df
# #   for (col_index in indices) {
# #     formatted_df[, col_index] <- color_tile_mean()(formatted_df[, col_index])
# #   }
# #   formattable(formatted_df)
# #   }
# #   else{
# #   p <- formattable(df, align=c("c", "c", "c"),list(
# #     col.names = as.character(col_names_bold)))
# #   }
# # }
# # c <- format_columns_by_indices(df)
# # length(c())
# # p <- formattable(df, align=c("c", "c", "c"),list(
# #   1=color_tile_mean(), 2=color_tile_mean(), col.names = as.character(col_names_bold)))


# #make_nice_table(df=df, rg_col_1 = 'microarray', outname = '/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/final_figures_and_tables/RT-qPCR_microarray_value_comparison.png')


# #    bad <- read.delim('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp.ngd')
# # table <- read.delim('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp_edited_no_symbols.ngd')
# # able <- table[,c(1,4,5,6)]
# # able <- rename(able, SEQ_ID = V2) 
# # na_table <- table[is.na(table$START),]
# # target_table <- drop_na(table)
# # target_table <- target_table %>% distinct(V2, .keep_all = TRUE)

# # lookup <- df2 %>% distinct(x, z)
# # merged <- merge(table, target_table, by = 'V2',all.x = TRUE)
# # write.table(merged, '~/Desktop/merged_ngd_tables.txt', quote = FALSE,row.names = FALSE,sep = '\t')

# # ngd.file <-read.delim('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp_andrea.ngd')



# # mock <- dir('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/mock/')
# # dex <- dir('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/dex/')
# # N <- length(mock)


# # for (i in 1:N) {
# #    if (i == 1) {pairs <- readPairs(mock[i],dex[i])
    
# #       }
# #    else {
# #      pairs <- readPairs(mock[i],dex[i],pairs)
# #       }
# #   }

# # pairs <- readPairs('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/mock/GSM2436560_561345A07_140403_450_532.pair','~/Documents/Ito_35S_microarray_data/GSE92729_RAW/dex/GSM2436563_561345A08_140403_450_532.pair')


# # ?file.path
# # ndf.file <-file.path('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp.ndf')
# # ngd.file <-file.path('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp.ngd')
# # pairs <- readDesign(ndf.file,ngd.file,pairs)

# ###Trichome_gene_anno 
# df <- read.csv('/Users/josephbeegan/Salba_RNA/genelists/trichome_GRN_salba.csv') %>% dplyr::select(gene_id, tair, gene_name)

# new_df <- df %>%
#   group_by(tair) %>%
#   summarise(gene_id = paste(gene_id, collapse = ","))

# new_df <- new_df %>% left_join(df, by = 'tair') %>% dplyr::select(-gene_id.y) %>% 
# dplyr::rename(TAIR_ID=tair, Salba_ID = gene_id.x, 'Symbol'=gene_name) %>% 
# distinct()

# make_nice_table(new_df,outname = '~/Salba_RNA/genelists/trichome_GRN_table.pdf')




# ###Making table of TCL1 novel gene counts
# samps <- s2c$Reads
# samps <- sort(samps)
# fpkm <- c(2.69,8.277,4.3,8.23,5.38,7.64,3.57,7.54,13.31,0,0,0,1.96,0,1.14,4.31,2.71,4.89,0.76,1.36,3.56,2.68,2.68,5.45,9.22,14.17,13.49,
# 6.60,7.70,7.11,19.94,22.42,22.87)



# g <- data_frame('Samples'=sort(samps), 'FPKM'=fpkm)

# g$exp <- substr(as.character(g$Samples), 1, nchar(as.character(g$Samples)) - 9)
# averages <- g %>%
#   group_by(exp) %>%
#   summarize(Avg_FPKM = mean(FPKM), 
#             StdDev = sd(FPKM), 
#             SE = StdDev / sqrt(n())) %>% mutate(exp=gsub('Early','Stage_5-7',exp))%>% mutate(exp=gsub('Mid','Stage_8-9',exp))%>% 
#             mutate(exp=gsub('Late','Stage_10-11',exp)) %>% arrange(exp)

# nb.cols <- 12
# mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)



# p<-ggplot(averages, aes(x = exp, y = Avg_FPKM, ymin = Avg_FPKM - SE, ymax = Avg_FPKM + SE, fill=exp)) +
#   geom_bar(stat = "identity",color="black", position=position_dodge()) +
#   geom_errorbar(width = 0.2, position = position_dodge(0.9)) +
#   labs(title = "Average FPKM",
#        x = "Tissue and Stage",
#        y = "Average FPKM")+
#        theme+
#   # scale_fill_brewer(palette = mycolors)+
#   scale_fill_manual(values = mycolors)+
#   # facet_wrap(~face)+
#   theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90),axis.line = element_line(colour = "black"))+
#    scale_x_discrete(limits = c(
#     "BUD-Stage_5-7", "BUD-Stage_8-9", "BUD-Stage_10-11",
#     "GYN-Stage_8-9", "GYN-Stage_10-11",
#     "STA-Stage_8-9", "STA-Stage_10-11",
#     "PET-Stage_8-9", "PET-Stage_10-11",
#     "SEP-Stage_8-9", "SEP-Stage_10-11"))+
#     guides(fill=guide_legend(title="Tissue and Stage"))+
#   ggtitle(paste("Average FPKM of Replicates", " Sal06g30330L"))

# out <- paste0('./gffcompare/basic_gffcompare/TCL1_FPKM_pattern.pdf')  
# print(out)
# ggsave(out,p)
