library(dplyr)
library(topGO)
library(KEGGREST)
library(org.At.tair.db)
library(GO.db)
library(ggplot2)
library(ggrepel)
library(wesanderson)
library(viridis)

194	Complete BUSCOs (C)			   
157	Complete and single-copy BUSCOs (S)	   
37	Complete and duplicated BUSCOs (D)	   
34	Fragmented BUSCOs (F)			   
27	Missing BUSCOs (M)			   
255	Total BUSCO groups searched	


Complete_BUSCOs= c(252, 247, 194)
Complete_and_single_copy_BUSCOs = c(96, 122, 157)
Complete_and_duplicated_BUSCOs=c(156,125,37)
Fragmented_BUSCOs=c(3,7,34)
Missing_BUSCOs=c(0,1,27)
mylist <- list(Complete_BUSCOs= c(252, 247, 173), Complete_and_single_copy_BUSCOs = c(96, 122, 117), Complete_and_duplicated_BUSCOs=c(156,125,56), Fragmented_BUSCOs=c(3,7,66), Missing_BUSCOs=c(0,1,16),Transcriptome = c('Trinity','SPAdesRNA','SOAP_aligned')) 
df <- as.data.frame(mylist)
df 
tdf <- transform.data.frame(df)
tdf
trinity_df <- data.frame(Title=c('Complete BUSCOs Trinity', "Complete and single-copy BUSCOs Trinity", "Complete and duplicated BUSCOs Trinty",'Fragmented BUSCOs Trinity','Missing BUSCOs Trinity','Complete BUSCOs SPAdesRNA', "Complete and single-copy BUSCOs SPAdesRNA", "Complete and duplicated BUSCOs SPAdesRNA",'Fragmented BUSCOs SPAdesRNA','Missing BUSCOs SPAdesRNA','Complete BUSCOs SOAP_denovo', "Complete and single-copy BUSCOs SOAP_denovo", "Complete and duplicated BUSCOs SOAP_denovo",'Fragmented BUSCOs SOAP_denovo','Missing BUSCOs SOAP_denovo'),
                         Num=c(252, 96, 156,3,0,247,122, 125,7,1,173,117,56,66,16), Transcriptome = c('Trinity','Trinity','Trinity','Trinity','Trinity','SPAdesRNA','SPAdesRNA','SPAdesRNA','SPAdesRNA','SPAdesRNA','SOAPdenovo','SOAPdenovo','SOAPdenovo','SOAPdenovo','SOAPdenovo'))
trinity_df
trinity_df <- trinity_df[order(trinity_df$Title),]
trinity_df
'SPAdesRNA_Complete BUSCOs', "SPAdesRNA_Complete and single-copy BUSCOs", "SPAdesRNA_Complete and duplicated BUSCOs",'SPAdesRNA_Fragmented BUSCOs','SPAdesRNA_Missing BUSCOs'),
len=c(247, 122, 125,7,1))
soapdenovo_aligned_df <- data.frame(dose=c('SOAP_aligned_Complete BUSCOs', "SOAP_aligned_Complete and single-copy BUSCOs", "SOAP_aligned_Complete and duplicated BUSCOs",'SOAP_aligned_Fragmented BUSCOs','SOAP_aligned_Missing BUSCOs'),
                                    len=c(173, 117, 56,66,16))




trinity_plot<-ggplot(data=trinity_df,aes(y=Title,x=Num,fill = Transcriptome)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 90))+
  geom_text(aes(label=Num),hjust = 2,color='black',size=5)
trinity_plot
spades_plot <- ggplot(data=spadesrna_df,aes(x=dose,y=len)) +
  geom_bar(stat="identity", fill = 'lightblue') +
  geom_text(aes(label=len), vjust=1.6, color="white", size=7.5)+
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 45)) +
  theme_minimal()
soapdenovo_aligned_plot <- ggplot(data=soapdenovo_aligned_df,aes(x=dose,y=len)) +
  geom_bar(stat="identity", fill = 'lightblue') +
  geom_text(aes(label=len), vjust=1.6, color="white", size=7.5)+
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 90)) +
  scale_colour_manual(values=c("Complete BUSCOs"="red", "Complete and single-copy BUSCOs"= "blue", "Complete and duplicated BUSCOs"="black", "Fragmented BUSCOs"="yellow",'Missing BUSCOs'= 'black'))
soapdenovo_aligned_plot
trinity_plot
spades_plot