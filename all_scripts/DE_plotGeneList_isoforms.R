DE_plotGeneList_isoforms <- function(counts, # dataframe with counts in long format + factors
                                     info, # info about the gene
                                     x,
                                     counts_col = "TPM",
                                     fill_boxplot,
                                     colour_jitter,
                                     ens_gene = "ens_gene",
                                     ext_gene = "ext_gene",
                                     description = "description",
                                     ann_db = NULL,
                                     save = "",
                                     dev = "pdf",
                                     th = theme()) {
  for (i in seq_len(nrow(info))) {
    ens <- info[i, ens_gene]
    ext <- info[i, ext_gene]
    desc <- info[i, description]
    tpm_factor <- dplyr::filter(counts, counts[ens_gene] == ens)

    p <- ggplot2::ggplot(tpm_factor, aes_string(x, counts_col)) +
      geom_boxplot(aes_string(fill = fill_boxplot), alpha = 0.3, outlier.shape = NA) +
      geom_jitter(aes_string(colour = colour_jitter),
        position = position_jitterdodge(
          jitter.width = 0.3,
          dodge.width = 0.9
        ),
        size = 2
      ) +
      th +
      theme(
        axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, colour = "black"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      ) +
      ggtitle(
        label = paste0("Ensembl ID = ", ens, " [ gene name = ", ext, " ]"),
        subtitle = desc
      )

    wh <- GenomicFeatures::genes(transcript_ann)[ens]
    wh <- range(wh, ignore.strand = T)
    q <- ggplot() +
      ggbio::geom_alignment(
        transcript_ann,
        which = wh,
        cds.rect.h = 0.10,
        arrow.rate = 0.045,
        label.color = "black"
      ) +
      th +
      scale_y_continuous(expand = c(0.2, 0)) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      )

    all <- ggpubr::ggarrange(p, q, ncol = 1, nrow = 2, heights = c(3, 1))
    ggsave(
      plot = all,
      file = paste0(save, ens, "_", ext, ".", dev),
      device = dev,
      width = 12,
      height = 9
    )

    print(paste("Saved:", paste0(save, ens, "_", ext, ".pdf")))
  }
}
