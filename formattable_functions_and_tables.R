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

#YOU NEED GOOGLE CHROME INSTALLED FOR THIS
#Function to export the table as png
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2){
      w <- as.htmlwidget(f, width = width, height = height)
      path <- html_print(w, background = background, viewer = NULL)
      url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
      webshot2::webshot(url,
              file = file,
              selector = ".formattable_widget",
              delay = delay)
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
          `background-color` = ifelse(x < 0 , '#f45a42', '#7eec78')) # Remember to change the colors!
  })}

#These two lines will make the titles of each column bold
col_names_bold <- sprintf("<b>%s</b>", df)
format_columns_by_indices <- function(df, indices=c()) {
  if(length(indices) >0){
  formatted_df <- df
  for (col_index in indices) {
    formatted_df[, col_index] <- color_tile_mean()(formatted_df[, col_index])
  }
  formattable(formatted_df, list(col_names_bold))
  }
  else{
  p <- formattable(df, align=c("c", "c", "c"),list(
    col.names = as.character(col_names_bold)))
  }
}





make_nice_table <- function(df,colnums=c(),outname){
  p <- format_columns_by_indices(df, indices = colnums)
  export_formattable(p, outname)
}

#make_nice_table(df=df, colnums = c(2,3), outname ='/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/final_figures_and_tables/RT-qPCR_microarray_value_comparison.png')



# df <- read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/final_figures_and_tables/RT-qPCR_microarray_value_comparison.csv')


# format_columns_by_indices <- function(df, indices=c()) {
#   if(length(indices) >0){
#   formatted_df <- df
#   for (col_index in indices) {
#     formatted_df[, col_index] <- color_tile_mean()(formatted_df[, col_index])
#   }
#   formattable(formatted_df)
#   }
#   else{
#   p <- formattable(df, align=c("c", "c", "c"),list(
#     col.names = as.character(col_names_bold)))
#   }
# }
# c <- format_columns_by_indices(df)
# length(c())
# p <- formattable(df, align=c("c", "c", "c"),list(
#   1=color_tile_mean(), 2=color_tile_mean(), col.names = as.character(col_names_bold)))


#make_nice_table(df=df, rg_col_1 = 'microarray', outname = '/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/SUP_microarray_results/Final_data/final_figures_and_tables/RT-qPCR_microarray_value_comparison.png')


#    bad <- read.delim('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp.ngd')
# table <- read.delim('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp_edited_no_symbols.ngd')
# able <- table[,c(1,4,5,6)]
# able <- rename(able, SEQ_ID = V2) 
# na_table <- table[is.na(table$START),]
# target_table <- drop_na(table)
# target_table <- target_table %>% distinct(V2, .keep_all = TRUE)

# lookup <- df2 %>% distinct(x, z)
# merged <- merge(table, target_table, by = 'V2',all.x = TRUE)
# write.table(merged, '~/Desktop/merged_ngd_tables.txt', quote = FALSE,row.names = FALSE,sep = '\t')

# ngd.file <-read.delim('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp_andrea.ngd')



# mock <- dir('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/mock/')
# dex <- dir('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/dex/')
# N <- length(mock)


# for (i in 1:N) {
#    if (i == 1) {pairs <- readPairs(mock[i],dex[i])
    
#       }
#    else {
#      pairs <- readPairs(mock[i],dex[i],pairs)
#       }
#   }

# pairs <- readPairs('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/mock/GSM2436560_561345A07_140403_450_532.pair','~/Documents/Ito_35S_microarray_data/GSE92729_RAW/dex/GSM2436563_561345A08_140403_450_532.pair')


# ?file.path
# ndf.file <-file.path('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp.ndf')
# ngd.file <-file.path('~/Documents/Ito_35S_microarray_data/GSE92729_RAW/GPL13970_100718_Athal_TAIR9_exp.ngd')
# pairs <- readDesign(ndf.file,ngd.file,pairs)
