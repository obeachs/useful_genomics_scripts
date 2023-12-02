library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(data.table)
library(BiocParallel)
library(parallel)
library(stringr)
library(tximport)
library(DRIMSeq)
library(GenomicRanges)

# get the number of threads of the machine
threads <- detectCores()
register(MulticoreParam(threads))

# set a cool color palette
tropical <- c("darkorange", "limegreen", "dodgerblue", "darkorchid")
palette(tropical)

# GGplot2 theme
th <- theme(
  axis.text.x = element_text(angle = 65, vjust = 0.5, colour = "black"),
  panel.background = element_blank(), panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)

# Add gene attributes from biomaRt
source("./Rlib/EnsemblAttributes.R")
t2g <- EnsemblAttributes(
  attributes = c(
    "ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name", "transcript_biotype",
    "description", "chromosome_name",
    "transcript_start", "transcript_end", "strand"
  )
) %>% dplyr::rename(
  feature_id = ensembl_transcript_id,
  gene_id = ensembl_gene_id,
  ext_gene = external_gene_name
)

# create a df with  directories for each sample and attributes
s2c <- read.table("./ids.csv",
  sep = ",",
  header = T,
  stringsAsFactors = T
) %>%
  filter(!str_detect(Sample, "FAIL")) %>%
  droplevels.data.frame() %>%
  mutate(path = paste0("./quantification_clean/quant/", Sample, "/quant.sf")) %>%
  mutate(DAI = factor(DAI))

# create a vector with  directories for each sample
kal_dirs <- s2c$path
names(kal_dirs) <- s2c$Sample

# find gene/transcript correlation to use in tximport.
tx2gene <- read.table("./assembly/ballgown/fl_0days_rep1_clean/t_data.ctab", sep = "\t", header = T, stringsAsFactors = FALSE) %>%
  dplyr::select(t_name, gene_id) %>%
  dplyr::rename(feature_id = t_name, gene_id = gene_id)

# Load a database of the transcript in the novel assembly
# I had to get rid of the few transcript with unknown strand position
transcript_ann <- GenomicFeatures::makeTxDbFromGFF("./assembly/merge_clean_fixstrand.gtf", format = "gtf")


# DIAGNOSTICS ________________________________________________________________________________________________

dir.create("diagnostics/gene_level/PCA/", recursive = T)
dir.create("diagnostics/transcript_level/PCA/", recursive = T)

source("./Rlib/DESeq2_plotstats.R")
source("./Rlib/PCA.R")

# gene level diagnostics
txi.gene <- tximport(kal_dirs,
  txOut = F,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "scaledTPM",
  varReduce = T
)


pca_plot(
  data = txi.gene$counts,
  factors = s2c,
  color = "DAI",
  th = th,
  text = T,
  save = "./diagnostics/gene_level/PCA/PCA_gene.pdf"
)

# trascript level diagnostics
txi.salmon <- tximport(
  kal_dirs,
  txOut = T,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "scaledTPM",
  varReduce = T
)

pca_plot(
  data = txi.salmon$counts,
  factors = s2c,
  color = "DAI",
  th = th,
  text = T,
  save = "./diagnostics/transcript_level/PCA/PCA_transcript.pdf"
)


# Contribution of each gene to PCs ___________________________________________________________________________

dir.create("diagnostics/gene_level/PC2_contributions/plots/", recursive = T)
dir.create("diagnostics/transcript_level/PC2_contributions/plots/", recursive = T)

source("./Rlib/PC_weights.R")
source("./Rlib/DE_plotGeneList_isoforms.R")

contribs_gene <- PC_top_contributors(txi_df = txi.gene$counts) %>%
  dplyr::rename("gene_id" = 1) %>%
  left_join(dplyr::select(t2g, gene_id, ext_gene, description)) %>%
  arrange(desc(PC2)) %>%
top_n(10,PC2) 

write.table(
  contribs_gene,
  "./diagnostics/gene_level/PC2_contributions/PC2_contributions_top10.csv",
  sep = ",",
  row.names = F
)

contribs_transcript <- PC_top_contributors(txi_df = txi.salmon$counts) %>%
  dplyr::rename("feature_id" = 1) %>%
  left_join(tx2gene) %>%
  left_join(dplyr::select(t2g, feature_id, gene_id, ext_gene, description)) %>%
  arrange(desc(PC2)) %>%
top_n(10,PC2) 

write.table(
  contribs_transcript,
  "./diagnostics/transcript_level/PC2_contributions/PC2_contributions_top10.csv",
  sep = ",",
  row.names = F
)

contribs_gene <- dplyr::select(contribs_gene, gene_id, ext_gene, description)
contribs_transcript <- unique(dplyr::select(contribs_transcript, gene_id, ext_gene, description))

plot_counts <- txi.salmon$counts %>%
  as.data.frame() %>%
  rownames_to_column("feature_id") %>%
  left_join(tx2gene) %>%
  pivot_longer(
    cols = c(-gene_id, -feature_id),
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  left_join(dplyr::select(s2c, -path))

contribs_gene_pc <- plot_counts %>%
  filter(gene_id %in% contribs_gene$gene_id)
contribs_transcript_pc <- plot_counts %>%
  filter(gene_id %in% contribs_transcript$gene_id)


DE_plotGeneList_isoforms(
  counts = contribs_gene_pc,
  info = contribs_gene,
  x = "feature_id",
  counts_col = "TPM",
  fill_boxplot = "DAI",
  colour_jitter = "DAI",
  ens_gene = "gene_id",
  ext_gene = "ext_gene",
  description = "description",
  ann_db = transcript_ann,
  save = "./diagnostics/gene_level/PC2_contributions/plots/",
  dev = "pdf",
  th = th
)

DE_plotGeneList_isoforms(
  counts = contribs_transcript_pc,
  info = contribs_transcript,
  x = "feature_id",
  counts_col = "TPM",
  fill_boxplot = "DAI",
  colour_jitter = "DAI",
  ens_gene = "gene_id",
  ext_gene = "ext_gene",
  description = "description",
  ann_db = transcript_ann,
  save = "./diagnostics/transcript_level/PC2_contributions/plots/",
  dev = "pdf",
  th = th
)

# PCA with 4days_rep3 removed ________________________________________________________________________________

kal_dirs <- kal_dirs[s2c$Sample != "4days_rep3"]
s2c <- filter(s2c, Sample != "4days_rep3") %>% droplevels.data.frame()

# gene level
txi.gene <- tximport(
  kal_dirs,
  txOut = F,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "scaledTPM",
  varReduce = T
)

pca_plot(
  data = txi.gene$counts,
  factors = s2c,
  color = "DAI",
  th = th,
  text = T,
  save = "./diagnostics/gene_level/PCA/PCA_gene_reduced.pdf"
)

# trascript level
txi.salmon <- tximport(
  kal_dirs,
  txOut = T,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "scaledTPM",
  varReduce = T
)

pca_plot(
  data = txi.salmon$counts,
  factors = s2c,
  color = "DAI",
  th = th,
  text = T,
  save = "./diagnostics/transcript_level/PCA/PCA_transcript_reduced.pdf"
)


# DRIMSeq analysis ___________________________________________________________________________________________

dir.create("DRIMSeq/plots/dtu_genes/", recursive = T)
dir.create("DRIMSeq/tables", recursive = T)

counts_drimseq <- txi.salmon$counts %>%
  as.data.frame() %>%
  rownames_to_column("feature_id") %>%
  left_join(tx2gene)
s2c <- dplyr::select(s2c, Sample, DAI) %>% dplyr::rename(sample_id = "Sample")
s2c$DAI <- as.character(s2c$DAI)

d <- dmDSdata(counts = counts_drimseq, samples = s2c)
plotData(d)

filtered <- dmFilter(d,
  min_samps_gene_expr = length(counts_drimseq) - 2,
  min_samps_feature_expr = length(counts_drimseq) - 2,
  min_samps_feature_prop = length(counts_drimseq) - 2,
  min_gene_expr = 50,
  min_feature_expr = 10,
  min_feature_prop = 0.1,
  run_gene_twice = FALSE
)


design_full <- model.matrix(~DAI, data = samples(d))

set.seed(123)
BPPARAM <- BiocParallel::MulticoreParam(workers = threads)

precision <- dmPrecision(d, design = design_full, BPPARAM = BPPARAM)
plotPrecision(precision)

fitted <- dmFit(precision,
  design = design_full, verbose = 1,
  BPPARAM = BPPARAM
)

design_null <- model.matrix(~1, data = samples(d))
tested <- dmTest(fitted, design = design_null, verbose = 1, BPPARAM = BPPARAM)

final_results <- tested@results_gene %>%
  filter(adj_pvalue < 0.05 & pvalue < 0.05) %>%
  left_join(dplyr::select(t2g, gene_id, ext_gene, transcript_biotype, description)) %>%
  unique.data.frame() %>%
  arrange(desc(adj_pvalue))

write.table(final_results, "./DRIMSeq/tables/DRIMseq_results_full.csv", sep = ",", row.names = F)


# Plot significant genes _____________________________________________________________________________________

fr_info <- final_results %>%
  dplyr::select(gene_id, ext_gene, description)
fr_pc <- plot_counts %>%
  filter(gene_id %in% final_results$gene_id)

source("./Rlib/DE_plotGeneList_isoforms.R")
DE_plotGeneList_isoforms(
  counts = fr_pc,
  info = fr_info,
  x = "feature_id",
  counts_col = "TPM",
  fill_boxplot = "DAI",
  colour_jitter = "DAI",
  ens_gene = "gene_id",
  ext_gene = "ext_gene",
  description = "description",
  ann_db = transcript_ann,
  save = "./DRIMSeq/plots/dtu_genes/",
  dev = "pdf",
  th = th
)
