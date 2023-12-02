get_significant_GO <- function(weighted_genelist,
                               gene_col = "gene_id",
                               weight_col = "pval",
                               ontology = "BP",
                               fun_selection = function(x) {
                                 return(x < 0.01)
                               },
                               save = "GO/") {
  gene_list <- weighted_genelist[, weight_col]
  names(gene_list) <- weighted_genelist[, gene_col]

  GOdata <- methods::new(
    "topGOdata",
    description = "BP analysis",
    ontology = ontology,
    allGenes = gene_list,
    nodeSize = 10,
    geneSelectionFun = fun_selection,
    annot = annFUN.org,
    mapping = "org.At.tair.db"
  )

  result <- topGO::runTest(
    GOdata,
    algorithm = "weight",
    statistic = "fisher"
  )

  tabf <- topGO::GenTable(
    GOdata,
    pvalue = result,
    topNodes = length(result@score),
    numChar = 120
  ) %>%
    dplyr::filter(pvalue < 0.05, Annotated > 10) %>%
    dplyr::arrange(pvalue)

  utils::write.csv(tabf, file = paste0(save, "GO_terms.csv"), row.names = F)

  # wantedNodes can contain the name of the nodes we want to highlight
  topGO::printGraph(GOdata, result, firstSigNodes = 10, useInfo = "all", fn.prefix = paste0(save, "GO"), pdfSW = T)
}
