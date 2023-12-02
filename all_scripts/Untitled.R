library(ggplot2)
library(formattable)
library(tidyr)
library(dplyr)
library(reactable)
library(kableExtra)


table <- read.delim('/Volumes/seagatedrive/downloads/hitdata_full_bnapus_clean2.txt_clean_text_finished')
table2 <- table[c(1,9),]




jelly_columns <- read.delim("/Volumes/seagatedrive/sinapis_assembly_shenanigans/pandas_jellyfish_25mer_reads_only_dataframe",sep = '\t')
plot(jelly_columns[1:100,], type="l")
ggplot(data=jelly_columns, aes(x=Reads)) + geom_histogram(binwidth=1000)
quast_SPAdes59 <- read.delim('/Volumes/seagatedrive/sinapis_assembly_shenanigans/SOAPdenvo59mer_run_ragtag_no_correction/ragtag_output/quast_results/latest/report.tsv')
quast_SOAP59 <- read.delim('/Volumes/seagatedrive/sinapis_assembly_shenanigans/SPAdes59mer_ragtag_BO_scaffold_no_read_corrections/ragtag_output/quast_results/latest/report.tsv')
quast_SOAP <- read.delim('/Volumes/seagatedrive/sinapis_assembly_shenanigans/original_SOAP_ragtag_augustus_and_blast/SOAP_ragtag_with_read_correction_quast_results/latest/report.tsv')
quast_SPAdes <- read.delim('/Volumes/seagatedrive/sinapis_assembly_shenanigans/SPAdes_sinapis_assembly_ragtag/ragtag_output/quast_results/latest/report.tsv')
doublequast <- full_join(quast_SOAP59,quast_SPAdes59)
tripequast <- full_join(quast_SOAP,doublequast)
fullquast <- full_join(quast_SPAdes,tripequast)
fullquast <- as.data.frame(fullquast)
fullquast %>%
  kbl() %>%
  kable_material(c("striped", "hover")) %>%
  row_spec(7, color = 'white',background = 'green') %>%
  row_spec(13, color = 'white',background = 'green') %>%
  row_spec(17, color = 'white',background = 'green')
fullquast2 <- fullquast[c(1, 7, 17, 21),]
fullquast2 %>%
  kbl() %>%
  kable_material(c("striped", "hover"))%>%
  kable_styling(font_size = 15)
formattable(doublequast, list(Assembly= color_text("orange",'white'),
SOAP59mer_scaffolds = color_tile('blue','green'), SPAdes59mer_scaffolds = color_tile('yellow','red')))
