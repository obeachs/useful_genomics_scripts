library(dplyr)
library(topGO)
library(tidry)
library(ggplot2)
library(ggrepel)
library(formattable)
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

#Read in the tables
table1 <- read.delim('/Volumes/seagatedrive/downloads/gl1_hitdata.txt')
tab#table2 
#table3 
#etc etc

#Reduce table sizes to only the data you want -  in this case just looking at the query name and the name of the 
#protein domain

table1 <- table1[,c(1,9)]
write.table(table1,'/Volumes/seagatedrive/downloads/gl1_hitdata_cleaned_no_format.txt',sep = '\t', quote = FALSE, row.names = FALSE)#where n is the column number
#do this for all tables
temp <- as.data.frame(as.matrix(table1))
formattable(temp_
ggplot(temp,aes(Short.name,Query)) +
  geom_tile(aes(fill=Short.name),color = 'black')
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)
#then you can simply merge the tables with tidyr and dyplr
table_full_merge <- dplyr::full_join(table1, table2, by = 'X') # where X is the name of the column that you want to merge
# by, in this case it would be the protein domain column name
table_inner_merge <- dplyr::full_join(table1, table2, by = 'X')
#This would be only the matching between the tables. 
#Repeat this with other tables, add the new tables to the ones you just made
#Once you have a table you're happy with you can make a cool table:
# https://rfortherestofus.com/2019/11/how-to-make-beautiful-tables-in-r/ this is a good resource.
#You can always check how the table looks by typing the variable name and hitting enter in the console below
#To check the columns we would use table1$columnname in the console (most things will be predicted by R)
