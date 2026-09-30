# Dot plot for GO, KEGG, or STRING enrichment results
# The count column may contain a number or text such as "20 of 8096".

plot_enrichment_dot <- function(
  data,
  description_col = "description",
  count_col = "count_in_network",
  fdr_col = "false_discovery_rate",
  top_n = 10,
  title = NULL,
  wrap_width = 35
) {
  required_columns <- c(description_col, count_col, fdr_col)
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  plot_data <- data.frame(
    description = as.character(data[[description_col]]),
    count_raw = data[[count_col]],
    fdr = suppressWarnings(as.numeric(data[[fdr_col]])),
    stringsAsFactors = FALSE
  )

  plot_data$count <- suppressWarnings(
    as.numeric(sub("^\\s*([0-9.]+).*$", "\\1", as.character(plot_data$count_raw)))
  )

  plot_data <- plot_data |>
    dplyr::filter(
      !is.na(description),
      nzchar(description),
      !is.na(count),
      !is.na(fdr)
    ) |>
    dplyr::arrange(fdr, dplyr::desc(count)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(
      description = stringr::str_wrap(description, width = wrap_width),
      description = stats::reorder(description, count)
    )

  if (nrow(plot_data) == 0) {
    stop("No valid rows remain after count and FDR conversion.")
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = count, y = description, size = count, color = fdr)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(
      low = "#e06663",
      high = "#3380bc"
    ) +
    ggplot2::labs(
      title = title,
      x = "Gene count",
      y = NULL,
      size = "Count",
      color = "FDR"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = "black",
        size = 18,
        face = "bold",
        hjust = 0.5
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14,
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        color = "black",
        size = 11,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(color = "black", size = 11),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

# Example
# install.packages(c("ggplot2", "dplyr", "stringr"))
# result_table <- read.csv(file.choose(), check.names = FALSE)
# p <- plot_enrichment_dot(
#   result_table,
#   description_col = "description",
#   count_col = "count_in_network",
#   fdr_col = "false_discovery_rate",
#   top_n = 10,
#   title = "GO Biological Process"
# )
# p
# ggplot2::ggsave("GO_enrichment_dotplot.png", p, width = 10, height = 7, dpi = 300)
