# GO and KEGG enrichment analysis for human gene symbols
# Required input: a character vector of HGNC gene symbols

run_go_kegg <- function(
  gene_symbols,
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.20,
  min_gs_size = 10,
  max_gs_size = 500
) {
  required_packages <- c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot",
    "ggplot2"
  )

  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Install required packages first: ",
      paste(missing_packages, collapse = ", ")
    )
  }

  gene_symbols <- unique(trimws(as.character(gene_symbols)))
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)]

  if (length(gene_symbols) == 0) {
    stop("No valid gene symbols were supplied.")
  }

  mapping <- clusterProfiler::bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )

  if (nrow(mapping) == 0) {
    stop("No gene symbols could be mapped to ENTREZID.")
  }

  entrez_ids <- unique(mapping$ENTREZID)

  go_result <- clusterProfiler::enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    readable = TRUE
  )

  kegg_result <- clusterProfiler::enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalue_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  )

  list(
    mapping = mapping,
    go = go_result,
    kegg = kegg_result
  )
}

plot_go <- function(go_result, show_category = 10) {
  result_table <- as.data.frame(go_result)

  if (nrow(result_table) == 0) {
    stop("The GO result contains no enriched terms.")
  }

  enrichplot::dotplot(
    go_result,
    split = "ONTOLOGY",
    showCategory = show_category
  ) +
    ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free") +
    ggplot2::scale_color_gradient(low = "#e06663", high = "#3380bc") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = NULL)
}

plot_kegg <- function(kegg_result, show_category = 10) {
  result_table <- as.data.frame(kegg_result)

  if (nrow(result_table) == 0) {
    stop("The KEGG result contains no enriched pathways.")
  }

  enrichplot::dotplot(
    kegg_result,
    showCategory = show_category
  ) +
    ggplot2::scale_color_gradient(low = "#e06663", high = "#3380bc") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = NULL)
}

# Example
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
# install.packages("ggplot2")
#
# deg <- read.csv(file.choose(), check.names = FALSE)
# enrichment <- run_go_kegg(deg$Gene.symbol)
# head(enrichment$mapping)
# plot_go(enrichment$go, show_category = 10)
# plot_kegg(enrichment$kegg, show_category = 10)
