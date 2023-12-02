library(dplyr)
library(biomaRt)
library(tidyr)
library(ggplot2)
library(tximport)
library(rhdf5)
library(tibble)
library(data.table)
library(BiocParallel)
library(parallel)
library(DESeq2)
library(vsn)
library(stringr)
library(topGO)
library(org.At.tair.db)
library(pheatmap)

# get the number of threads of the machine
t <- detectCores()

register(MulticoreParam(t))

# set a cool color palette
tropical <- c("darkorange", "limegreen", "dodgerblue", "darkorchid")
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

# create a file with  directories for each sample and informations
s2c <- read.table("./ids.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = TRUE
) %>%
  mutate(path = paste0("./ballgown/", Reads, "/t_data.ctab")) %>%
  mutate(Sample = paste0(Genotype, "_", Treatment, "_", DAI, "DAI", "_Rep", Replica)) %>%
  mutate(Sample_pca = paste0(Genotype, "_Rep", Replica)) %>%
  mutate(DAI = as.factor(DAI), Replica = as.factor(Replica))

# Add gene attributes from biomaRt
source("./DE_scripts/EnsemblAttributes.R")
t2g <- EnsemblAttributes(
  attributes = c("transcript_biotype", "description", "chromosome_name")
)

# Import files from Stringtie output into R using tximport
file <- file.path(s2c$path[1])
tx2gene <- read.table(file, header = TRUE) %>%
  dplyr::select(t_name, gene_id)

txi <- tximport(
  s2c$path,
  type = "stringtie",
  txOut = FALSE,
  tx2gene = tx2gene,
  countsFromAbundance = "no"
)

# AG
s2c_AG <- s2c %>%
  filter(Genotype == "AG") %>%
  droplevels.data.frame()
txi_AG <- tximport(
  s2c_AG$path,
  type = "stringtie",
  txOut = FALSE,
  tx2gene = tx2gene,
  countsFromAbundance = "no"
)


# AP3
s2c_AP3 <- s2c %>%
  filter(Genotype == "AP3") %>%
  droplevels.data.frame()
txi_AP3 <- tximport(
  s2c_AP3$path,
  type = "stringtie",
  txOut = FALSE,
  tx2gene = tx2gene,
  countsFromAbundance = "no"
)


# Set seed for reproducibility -------------------------------------------
set.seed(33)

# Full dataset -----------------------------------------------------------
dir.create("plots/full", recursive = TRUE)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = s2c,
  design = ~ Genotype + Treatment + DAI
)

# DESeq2 analysis
dds <- DESeq(dds)

# Plot the dataset (healthcheck)
source("./DE_scripts/DESeq2_plotstats.R")
DESeq2_plotstats(
  dds,
  condition = c("Genotype"),
  sample_column = "Sample",
  save = "plots/full/",
  dev = "png",
  th = theme
)

# Full dataset: better PCA -----------------------------------------------
source("./DE_scripts/PCA.R")
tmp <- assay(dds)
colnames(tmp) <- s2c$Sample
pca_plot(
  data = tmp,
  factors = s2c,
  color = "DAI",
  shape = "Treatment",
  sample_names = "Sample_pca",
  th = theme,
  text = TRUE,
  save = "./plots/full/ag-ap3_full_PCA.pdf"
)

# AG dataset -----------------------------------------------------------------
dir.create("plots/AG", recursive = TRUE)

dds_AG <- DESeqDataSetFromTximport(
  txi_AG,
  colData = s2c_AG,
  design = ~ Treatment + DAI
)

# DESeq2 analysis
dds_AG <- DESeq(dds_AG, test = "LRT", reduced = ~DAI)

# Plot the dataset (healthcheck)
DESeq2_plotstats(
  dds_AG,
  condition = c("Treatment", "DAI"),
  sample_column = "Sample",
  save = "plots/AG/",
  dev = "png",
  th = theme
)

# Better PCA -------------------------------------------------------------
tmp <- assay(dds_AG)
colnames(tmp) <- s2c_AG$Sample
pca_plot(
  data = tmp,
  factors = s2c_AG,
  color = "DAI",
  shape = "Treatment",
  sample_names = "Sample_pca",
  th = theme,
  text = TRUE,
  save = "./plots/AG/ag-ap3_ag_PCA.pdf"
)

# AP3 dataset -----------------------------------------------------------------
dir.create("plots/AP3", recursive = TRUE)

dds_AP3 <- DESeqDataSetFromTximport(
  txi_AP3,
  colData = s2c_AP3,
  design = ~ Treatment + DAI
)

# DESeq2 analysis
dds_AP3 <- DESeq(dds_AP3, test = "LRT", reduced = ~DAI)

# Plot the dataset (healthcheck)
DESeq2_plotstats(
  dds_AP3,
  condition = c("Treatment", "DAI"),
  sample_column = "Sample",
  save = "plots/AP3/",
  dev = "png",
  th = theme
)

# save.image(file = "RData_afterddscreation.obj")
# load("RData_afterddscreation.obj")

# Better PCA -------------------------------------------------------------
tmp <- assay(dds_AP3)
colnames(tmp) <- s2c_AP3$Sample
pca_plot(
  data = tmp,
  factors = s2c_AP3,
  color = "DAI",
  shape = "Treatment",
  sample_names = "Sample_pca",
  th = theme,
  text = TRUE,
  save = "./plots/AP3/ag-ap3_ap3_PCA.pdf"
)

# Generate results ------------------------------------------------------------
source("./DE_scripts/DESeq2_saveresults.R")

dir.create("DE_results/AG/", recursive = TRUE)
dir.create("DE_results/AP3/", recursive = TRUE)

DESeq2_saveresults(
  dds_AG,
  condition = "Treatment",
  control = "Mock",
  info = dplyr::rename(t2g, gene_id = ens_gene),
  save = "DE_results/AG/",
  pval = 0.01,
  qval = 0.01,
  FC = 0,
  shrink = TRUE
)

DESeq2_saveresults(
  dds_AP3,
  condition = "Treatment",
  control = "Mock",
  info = dplyr::rename(t2g, gene_id = ens_gene),
  save = "DE_results/AP3/",
  pval = 0.01,
  qval = 0.01,
  FC = 0,
  shrink = TRUE
)

masterfile <- assay(dds) %>%
  as.data.frame()
colnames(masterfile) <- s2c$Sample
masterfile <- rownames_to_column(masterfile, "gene_id")


# Compare DEX over Mock for each timepoint and genotype ----------------------
for (genotype in c("AG", "AP3")) {
  for (day in c("3", "5", "8")) {
    #dir.create(paste0("plots/", genotype, "/", day, "DAI"), recursive = TRUE)
    s2c_tmp <- s2c %>%
      filter(Genotype == genotype) %>%
      filter(DAI == day) %>%
      droplevels.data.frame()

    txi_tmp <- tximport(
      s2c_tmp$path,
      type = "stringtie",
      txOut = FALSE,
      tx2gene = tx2gene,
      countsFromAbundance = "no"
    )

    dds_tmp <- DESeqDataSetFromTximport(
      txi_tmp,
      colData = s2c_tmp,
      design = ~Treatment
    )

    # DESeq2 analysis
    dds_tmp <- DESeq(dds_tmp, test = "LRT", reduced = ~1)


    res <- DESeq2::results(dds_tmp, contrast = c("Treatment", "DEX", "Mock"))
    res <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, padj, log2FoldChange, pvalue) %>%
    dplyr::rename(!!paste0("padj_", genotype, "_", day) := padj, 
        !!paste0("log2FoldChange", genotype, "_", day) := log2FoldChange,
        !!paste0("pvalue", genotype, "_", day) := pvalue)

    masterfile <- left_join(masterfile, res)

    ## Plot the dataset (healthcheck)
    #DESeq2_plotstats(
    #  dds_tmp,
    #  condition = c("Treatment"),
    #  sample_column = "Sample",
    #  save = paste0("plots/", genotype, "/", day, "DAI/"),
    #  dev = "png",
    #  th = theme
    #)

    #dir.create(paste0("DE_results/", genotype, "/", day, "DAI"), recursive = TRUE)

    #DESeq2_saveresults(
    #  dds_tmp,
    #  condition = "Treatment",
    #  control = "Mock",
    #  info = dplyr::rename(t2g, gene_id = ens_gene),
    #  save = paste0("DE_results/", genotype, "/", day, "DAI/"),
    #  pval = 0.01,
    #  qval = 0.01,
    #  FC = 0,
    #  shrink = TRUE
    #)
  }
}

masterfile  <- left_join(masterfile, t2g, by = c("gene_id" = "ens_gene"))
#write.csv(masterfile, file = "./masterfile_AG-AP3.csv", quote = FALSE, row.names = FALSE)
write.csv(masterfile, file = "./masterfile_AG-AP3_pval.csv", quote = FALSE, row.names = FALSE)

# Group all the results per DAI in one big dataset ---------------------------
res <- list.dirs(recursive = TRUE) %>%
  str_subset(pattern = "DAI/") %>%
  str_subset(pattern = "GO", negate = TRUE)

dat <- data.frame()
for (dir in res) {
  gene <- stringr::str_match(string = dir, pattern = "DE_results\\/(.*)\\/\\dDAI")[2]
  dai <- stringr::str_extract(string = dir, pattern = "\\dDAI")
  tmp <- read.table(paste0(dir, "/DEX.csv"), sep = ",", header = TRUE, stringsAsFactors = TRUE)
  tmp <- tibble::add_column(tmp, DAI_timepoint = dai, .before = 2)
  tmp <- tibble::add_column(tmp, Gene_knockdown = gene, .before = 2)
  dat <- rbind.data.frame(dat, tmp)
}
dat <- dplyr::arrange(dat, gene_id)

root <- "./DE_results/full_datasets_tables/"
dir.create(root, recursive = TRUE)
write.csv(dat, file = paste0(root, "AG-AP3_allDAI.csv"), row.names = FALSE)
write.csv(dplyr::filter(dat, Gene_knockdown = "AG"), file = paste0(root, "AG_allDAI.csv"), row.names = FALSE)
write.csv(dplyr::filter(dat, Gene_knockdown = "AP3"), file = paste0(root, "AP3_allDAI.csv"), row.names = FALSE)


# Plot key genes -------------------------------------------------------------
source("./DE_scripts/DE_plotGeneList.R")
dir.create("./plots/test_genes/", recursive = TRUE)

gene_list <- read.table("./ext_data/test_genes.csv", sep = ",", header = TRUE)$ens_gene

# HAS to have "Treatment" column (for now).
DE_plotGeneList(
  dds,
  t2g,
  gene_list = gene_list,
  condition = "Genotype",
  treatment = "Treatment",
  series = "DAI",
  save = "./plots/test_genes/",
  th = theme
)

dir.create("./plots/AG/lrt_Treatment/", recursive = TRUE)

gene_list <- read.table("./DE_results/AG/reference_MOCK/DEX.csv", sep = ",", header = TRUE)$gene_id

# HAS to have "Treatment" column (for now).
DE_plotGeneList(
  dds,
  t2g,
  gene_list = gene_list,
  condition = "Genotype",
  treatment = "Treatment",
  series = "DAI",
  save = "./plots/AG/lrt_Treatment/",
  th = theme
)

dir.create("./plots/AP3/lrt_Treatment/", recursive = TRUE)

gene_list <- read.table("./DE_results/AP3/reference_MOCK/DEX.csv", sep = ",", header = TRUE)$gene_id

# HAS to have "Treatment" column (for now).
DE_plotGeneList(
  dds,
  t2g,
  gene_list = gene_list,
  condition = "Genotype",
  treatment = "Treatment",
  series = "DAI",
  save = "./plots/AP3/lrt_Treatment/",
  th = theme
)
