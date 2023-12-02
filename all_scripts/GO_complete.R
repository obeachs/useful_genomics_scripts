library(topGO)
library(dplyr)
library(ggplot2)
source('/Volumes/sesame/ALP_Omics/ChIP/scripts/annotate_peaks.R')
#KG samples 


files <- list.files('~/Downloads/annotatedchipwithgaplesscalling(1)/kg_transfer/annnotated/', full.names = T, pattern = 'annotated')
for (file in files){
  p <- read.csv(file) %>% dplyr::filter(FDR_threshold < 0.01)
  outname <- basename(tools::file_path_sans_ext(file))
  write.csv(paste0('~/Downloads/annotatedchipwithgaplesscalling(1)/kg_transfer/', outname, 'FDR_0.01.csv'),quote=F, row.names=F)
}

for (file in files){
  location <- unlist(gregexpr('IP', basename(file)))[1]
  print(file)
  outname <-substr(basename(file),1,location + 1)
  outname <- paste(outname, '_annotated_nearest_option',sep = '')
  fuck <- read.csv(file, sep = '\t')
  print(names(fuck))
  annotate_peaks(file,out_file_prefix = paste('~/Downloads/annotatedchipwithgaplesscalling(1)/kg_transfer/',outname, sep = ''),
  chr_col = 'chrom', start_col = 'start', end_col = 'end', sep = ',')
}

# for (file in files){
#   gene <- basename(tools::file_path_sans_ext(file))
#   gene <- gsub('annotated','',gene)
#   peaks <- read.csv(file) %>% dplyr::select(feature)
#   peaks$pval <- 0.01
#   print(length(peaks$feature))
#   #get_significant_GO(peaks, gene_col = 'feature', save = paste('~/Downloads/annotatedchipwithgaplesscalling(1)/',name))
# }

test <- read.csv('~/Downloads/annotatedchipwithgaplesscalling(1)//kg_transfer/annnotated/AG_3d_GOGO_terms.csv')
View(test)
#Input should be list of genes and respective pvalues
#'save' option is just a directory of where to save the results
get_significant_GO <- function(weighted_genelist,
                               gene_col = "gene_id",
                               weight_col = "pval",
                               ontology = "BP",
                               fun_selection = function(x) {
                                 return(x < 0.01)
                               },
                               save = "GO/") {
  
  #formatting the data in a way that is readable for topGO
  #Every column is a gene name and the row of that column is its pvalue
  gene_list <- weighted_genelist[, weight_col]
  names(gene_list) <- weighted_genelist[, gene_col]
  gene_list[is.na(gene_list)] <- 1

  #Creating the topGO object using TAIR
  GOdata <- methods::new(
    "topGOdata",
    description = "BP analysis",
    ontology = ontology,
    allGenes = gene_list,
    nodeSize = 5,
    geneSelectionFun = fun_selection,
    annot = annFUN.org,
    mapping = "org.At.tair.db"
  )

#Getting significantly enriched GO terms with Fisher test
  resultFI <- topGO::runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  # resultGL <- topGO::runTest(
  #  GOdata,
  #  algorithm = "weight01",
  #  statistic = "globaltest"
  # )

#Generating a table of the GO term results and writing it to CSV files 
#Taking only significant values
 tabf <- topGO::GenTable(
    GOdata,
    topNodes = max(length(resultFI@score)),
    numChar = 120,
    Fisher = resultFI,
    orderBy = "Fisher",
    ranksOf = "Fisher"
  ) 
  print(head(tabf))
  tabf <- dplyr::filter(tabf,Fisher < 0.05) %>%
    dplyr::filter(Significant > 1) %>%
    dplyr::arrange(Fisher)

    utils::write.csv(tabf, file = paste0(save, "GO_terms.csv"), row.names = F)
    if (length(tabf$Fisher) > 100){
      tabf_short <- dplyr::arrange(tabf, desc(Annotated)) %>% slice(1:100)
      utils::write.csv(tabf_short, file = paste0(save, "GO_terms_top100.csv"), row.names = F)
    }

  # wantedNodes can contain the name of the nodes we want to highlight
  save(GOdata, file = paste0(save, "GOdata.obj"))
   dir.create(paste0(save, "genes_in_GO/"), recursive = TRUE)


  #Making a file showing the GO_terms and the genes associated
  genes <- dplyr::filter(weighted_genelist, FDR_threshold < 0.01)$feature
  GO2genes <- topGO::genesInTerm(GOdata)
  GO2genes_subset <- lapply(GO2genes, function(x) x[x %in% genes])
  GO2genes_subset <- GO2genes_subset[lapply(GO2genes_subset, length) > 0]
  tabf <- read.csv(paste0(save, "GO_terms.csv")) %>%
    dplyr::mutate(Term = stringr::str_replace_all(Term, "\\s+", "_")) %>%
    dplyr::mutate(Term = stringr::str_replace_all(Term, "\\/", "-")) %>%
    dplyr::mutate(Term = substr(Term, 1, 40))
  print(head(genes))
  for (go in tabf$GO.ID) {
    go_descr <- dplyr::filter(tabf, GO.ID == go)$Term
    print(go_descr)
    genes <- data.frame(gene_id = GO2genes_subset[[go]]) %>%
      utils::write.csv(file = paste0(save, "genes_in_GO/", go, "_", go_descr, ".csv"), row.names = FALSE)
    print(genes)
    print(paste0(save, "genes_in_GO/", go, "_", go_descr, ".csv"))
  #This generates the ugly plots/flowcharts
  #topGO::printGraph(GOdata, resultFI, firstSigNodes = 5, useInfo = "all", fn.prefix = paste0(save, "FI_flow_chart_first5"), pdfSW = TRUE)
  }
}


#Need to understand how the kmeans clustering works.

for (file in files){
  location <- unlist(gregexpr('d_', basename(file)))[1]
  outname <-substr(basename(file),1,location)
  outname <- paste(dirname(file),'/',outname,sep = '')
  genelist <- read.csv(file) %>% dplyr::select(feature,FDR_threshold) %>% dplyr::filter(FDR_threshold < 0.01)
  get_significant_GO(genelist,gene_col = 'feature',weight_col = 'FDR_threshold', save = paste(outname,'GO',sep = '_'))
}


test <- read.csv('/Users/josephbeegan/Downloads/annotatedchipwithgaplesscalling(1)/kg_transfer/AG_3d_IP_annotated_nearest_option.csv')
length(test$FDR_threshold[test$FDR_threshold < 0.01])



# set a cool color palette
tropical <- c("darkorange", "limegreen", "dodgerblue", "darkorchid", "lightcoral")
palette(tropical)

# GGplot2 theme
theme <- theme(
  axis.text.x = element_text(colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)


# All significant GOs ----------------------------------------------------------

dir.create("~/Downloads/annotatedchipwithgaplesscalling(1)/kg_transfer/annnotated/all_GO", recursive = TRUE)

gos <- data.frame(GO.ID = c(), Term = c())
for (gene in c("AG", "AP3")) {
  for (dai in c("3d", "5d", "8d")) {
    dat_tmp <- read.csv(paste0("~/Downloads/annotatedchipwithgaplesscalling(1)/kg_transfer/annnotated/top_100/", gene, "_", dai,'_GOGO_terms_top100.csv'), header = TRUE)
    gos <- rbind.data.frame(gos, dat_tmp[, 1:2])
  }
}
gos <- unique.data.frame(gos)

gos_df <- data.frame()
for (gene in c("AG", "AP3")) {
  for (dai in c("3d", "5d", "8d")) {
    dat_tmp <- read.csv(paste0("~/Downloads/annotatedchipwithgaplesscalling(1)/kg_transfer/annnotated/top_100/", gene, "_", dai,'_GOGO_terms_top100.csv'), header = TRUE)
    dat_tmp <- left_join(gos, dat_tmp) %>%
      tibble::add_column(Gene = gene, DAI = dai, .before = 1)
    gos_df <- rbind.data.frame(gos_df, dat_tmp)
  }
}

# Plot significant GOs as heatmaps ---------------------------------------------

# Fisher
gos_df_fisher <- gos_df %>%
  group_by(GO.ID) %>%
  filter(any(Fisher <= 0.001, na.rm = TRUE)) %>%
  ungroup()
gos_df_fisher$Fisher[is.na(gos_df_fisher$Fisher)] <- 1.0
gos_df_fisher$Fisher[gos_df_fisher$Fisher == "< 1e-30"] <- 0.0000001
gos_df_fisher$Fisher <- as.numeric(gos_df_fisher$Fisher)

gos_clust <- gos_df_fisher %>%
  dplyr::select(Gene, DAI, GO.ID, Fisher) %>%
  tidyr::pivot_wider(names_from = c(Gene, DAI), values_from = Fisher) %>%
  tibble::column_to_rownames("GO.ID")
tmp <- data.frame(kmeans(as.matrix(gos_clust), 5, 1000)$cluster) %>%
  tibble::rownames_to_column("GO.ID")
colnames(tmp) <- c("GO.ID", "cluster")


gos_df_fisher2 <- gos_df_fisher %>%
  left_join(tmp) %>%
  arrange(cluster, GO.ID) %>%
  mutate(cluster = factor(cluster)) 
gos_df_fisher2$Fisher[is.na(gos_df_fisher2$Fisher)] <- 1.0
gos_df_fisher2$Fisher <- as.numeric(gos_df_fisher2$Fisher)
gos_df_fisher2$Fisher[is.na(gos_df_fisher2$Fisher)] <- 1.0
gos_df_fisher2$rfisher <- log2(gos_df_fisher2$Fisher)





p <- ggplot(
  mutate(gos_df_fisher2, Term = paste0(stringr::str_sub(Term, end = 50))),
  aes(Gene, reorder(Term, as.numeric(cluster)))
) +
  facet_wrap(~Gene, ncol = 3) +
  geom_tile(aes(fill = cluster, alpha = rfisher ), color = "black", size = 0.7) +
  theme_minimal() +
  # scale_fill_gradient(high = "white", low = tropical[1]) +
  scale_fill_manual(values = palette()) +
  scale_alpha(range = c(1, 0)) +
  theme(
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_blank()
  )
ggsave(p, file = "~/Desktop/GO_terms_heatmap_Fisher.pdf", device = "pdf", scale = 0.7)
