#' Generate a Volcano Plot for DEG Analysis
#'
#' Creates a volcano plot from a differential gene expression result table.
#'
#' @param data A data frame containing gene, log2FoldChange, and pvalue columns.
#' @param log2fc_cutoff Numeric absolute log2 fold-change threshold.
#' @param pvalue_cutoff Numeric p-value threshold.
#' @param label_genes Logical; label significant genes when TRUE.
#' @param label_top_n Maximum number of significant genes to label, ordered by
#'   p-value.
#' @return A ggplot2 object.
#' @export
generate_volcano_plot <- function(
  data,
  log2fc_cutoff = 1,
  pvalue_cutoff = 0.05,
  label_genes = TRUE,
  label_top_n = 20
) {
  required_columns <- c("gene", "log2FoldChange", "pvalue")
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (!is.numeric(data$log2FoldChange) || !is.numeric(data$pvalue)) {
    stop("log2FoldChange and pvalue must be numeric columns.")
  }

  plot_data <- dplyr::mutate(
    data,
    pvalue_plot = pmax(pvalue, .Machine$double.xmin),
    neg_log10_p = -log10(pvalue_plot),
    Significance = dplyr::case_when(
      pvalue < pvalue_cutoff & log2FoldChange > log2fc_cutoff ~ "Up",
      pvalue < pvalue_cutoff & log2FoldChange < -log2fc_cutoff ~ "Down",
      TRUE ~ "NotSignif"
    )
  )

  significant_genes <- plot_data |>
    dplyr::filter(Significance != "NotSignif", !is.na(gene)) |>
    dplyr::arrange(pvalue) |>
    dplyr::slice_head(n = label_top_n)

  volcano_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = log2FoldChange, y = neg_log10_p, fill = Significance)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2, shape = 21, color = "black") +
    ggplot2::scale_fill_manual(
      values = c(Up = "red", Down = "blue", NotSignif = "grey")
    ) +
    ggplot2::geom_vline(
      xintercept = c(-log2fc_cutoff, log2fc_cutoff),
      color = "black",
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(pvalue_cutoff),
      color = "black",
      linetype = "dashed"
    ) +
    ggplot2::labs(
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      fill = "Significance"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text = ggplot2::element_text(face = "bold", size = 12),
      legend.title = ggplot2::element_text(face = "bold")
    )

  if (isTRUE(label_genes) && nrow(significant_genes) > 0) {
    volcano_plot <- volcano_plot +
      ggrepel::geom_label_repel(
        data = significant_genes,
        ggplot2::aes(label = gene),
        box.padding = 0.5,
        point.padding = 0.5,
        segment.color = "grey50",
        max.overlaps = Inf,
        fontface = "bold",
        size = 4,
        fill = "white"
      )
  }

  volcano_plot
}
