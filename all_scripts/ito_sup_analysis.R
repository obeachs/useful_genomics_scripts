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



'3.',
'0.2972383844',
'0.02491946464',
'1.558177299',
'1.316056163',
'-1.23908761',
'1.016307947',
'-0.2988956416',
'-3.777182896',
'4.408526002',
'0.7687614657',
'0.6975302468',
'-1.127986539',
'0.2014896141',
'0.2297677499',
'0.5603306287'


library(formattable)
color_formatter <- formattable::formatter(
  "span",
  style = x ~ style(
    color = 'white',
    'background-color' =
      ifelse(x >0, "green", "red")
            ))
formattable::formattable(df, list(
  RTqPCR = color_formatter,
  microarray=color_formatter,
   = formattable::color_tile("transparent","lightgreen")))


p <- formattable(GOI_df, align=c("c", "c", "c"),list(
  logFC=color_tile_mean()
  )
)
export_formattable(
  p,
  '/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GOI_DEGS',
  width = "100%",
  height = "100%",
  background = "white",
  delay = 10


df <- data.frame (Gene=c('PAP85', 'PMP','AT3G20280','PBP1','PMT5','CML1','CYP77A5P','PAP1-L','SUP','URO','AP3','FOA1','ESR2','WOX12','AHP6','AG'),
                  RTqPCR  = c('3.56828422',
'0.2972383844',
'0.02491946464',
'1.558177299',
'1.316056163',
'-1.23908761',
'1.016307947',
'-0.2988956416',
'-3.777182896',
'4.408526002',
'0.7687614657',
'0.6975302468',
'-1.127986539',
'0.2014896141',
'0.2297677499',
'0.5603306287'),
                  microarray = c("-2.01927808109384", "-1.31616728373971","-3.30096730095337","-1.07953621372366","1.16219358547157","-1.03898315618023",
                  "-1.1429516757837","-0.265338168666391","-3.49514743578471","1.14315146863367","-0.0872530992218845",
                  "0.0129064047038758","0.595337703542807","-1.44747860122867","0.327121921070924","0.0531503899772297")
                  )






deg_tab <- read.csv('/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/dex_v_mock_35S-SUP-GR_DEGs.tsv', sep = '\t')

genes <- deg_tab$ORF
GOI <- c('AT3G23130','AT3G23140','AT4G18960','AT3G54340')
tair_mart <- useMart(biomart='plants_mart',host = 'plants.ensembl.org', dataset='athaliana_eg_gene')
keytypes(org.At.tair.db)
mapIds(org.At.tair.db, keys = genes,
  column = c('SYMBOL'), keytype = 'TAIR')

eep <- biomaRt::select(org.At.tair.db, keys = genes,
  column = c('SYMBOL'), keytype = 'TAIR') %>% dplyr::rename(ORF=TAIR)


full_df <- left_join(deg_tab, eep,by='ORF') %>% dplyr::select(ORF, SYMBOL, logFC, adj.P.Val)
color_tile_mean <- function (...) {
  formatter("span", style = function(x) {
    style(display = "block",
          padding = "0 4px", 
          `border-radius` = "4px", 
          `background-color` = ifelse(x < 0 , '#f45a42', '#7eec78')) # Remember to change the colors!
  })}

GOI_df <- dplyr::filter(full_df, ORF %in% GOI)

p <- formattable(df, align=c("c", "c", "c"),list(
  RTqPCR=color_tile_mean(),
  microarray=color_tile_mean()
  )
)
export_formattable(
  p,
  '/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GOI_DEGS',
  width = "100%",
  height = "100%",
  background = "white",
  delay = 10
)

export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2){
      w <- as.htmlwidget(f, width = width, height = height)
      path <- html_print(w, background = background, viewer = NULL)
      url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
      webshot(url,
              file = file,
              selector = ".formattable_widget",
              delay = delay)
    }

export_formattable(p, '/Volumes/sesame/joerecovery/Project_folder/microarray_SUP/Ito_35S_microarray_data/GOI_DEGs.png')



customGreen0 = "#DeF7E9"

customGreen = "#71CA97"

customRed = "#ff7f7f"

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
