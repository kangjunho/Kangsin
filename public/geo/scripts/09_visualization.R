# Step 09. Prepare gene lists for GO and KEGG over-representation analysis
# Input files are created in Steps 07 and 08.

# 1. Install and load annotation packages ------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("AnnotationDbi", "org.Hs.eg.db")
missing_bioc <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

library(AnnotationDbi)
library(org.Hs.eg.db)

# 2. Load DEG and probe annotation results -----------------------------------
deg_dir <- "results/deg"
data_dir <- "data_processed"
output_dir <- "results/enrichment_input"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

deg <- read.csv(
  file.path(deg_dir, "GSE54568_limma_annotated.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
probe_map <- read.csv(
  file.path(data_dir, "probe_to_gene_valid.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_deg_cols <- c(
  "PROBE_ID", "logFC", "P.Value", "adj.P.Val",
  "gene_symbol", "status"
)
if (!all(required_deg_cols %in% colnames(deg))) {
  stop(
    "DEG 결과에 필요한 열이 없습니다: ",
    paste(setdiff(required_deg_cols, colnames(deg)), collapse = ", ")
  )
}

# 3. Clean symbols and prepare the measured background ------------------------
clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

deg$gene_symbol <- clean_symbol(deg$gene_symbol)
probe_map$gene_symbol <- clean_symbol(probe_map$gene_symbol)

# Background: all unique genes represented by analyzed, valid probes.
background_symbols <- sort(unique(stats::na.omit(
  probe_map$gene_symbol[
    probe_map$probe_id %in% deg$PROBE_ID
  ]
)))

# 4. Build direction-specific DEG symbol lists --------------------------------
up_symbols_raw <- sort(unique(stats::na.omit(
  deg$gene_symbol[deg$status == "Up"]
)))
down_symbols_raw <- sort(unique(stats::na.omit(
  deg$gene_symbol[deg$status == "Down"]
)))

# A gene may have different probes classified in opposite directions.
direction_conflicts <- intersect(up_symbols_raw, down_symbols_raw)

# Conservative beginner workflow: exclude conflicting symbols from both lists.
up_symbols <- setdiff(up_symbols_raw, direction_conflicts)
down_symbols <- setdiff(down_symbols_raw, direction_conflicts)
all_deg_symbols <- sort(unique(c(up_symbols, down_symbols)))

if (length(up_symbols) == 0) {
  warning("Up gene symbol이 없습니다. Up ORA를 강제로 수행하지 않습니다.")
}
if (length(down_symbols) == 0) {
  warning("Down gene symbol이 없습니다. Down ORA를 강제로 수행하지 않습니다.")
}

# 5. Convert SYMBOL to ENTREZID without discarding one-to-many mappings --------
symbols_to_entrez <- function(symbols, list_name) {
  symbols <- sort(unique(stats::na.omit(clean_symbol(symbols))))

  if (length(symbols) == 0) {
    return(data.frame(
      list = character(),
      SYMBOL = character(),
      ENTREZID = character(),
      stringsAsFactors = FALSE
    ))
  }

  converted <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = "ENTREZID"
  )
  converted <- converted[
    !is.na(converted$ENTREZID) &
      nzchar(converted$ENTREZID),
    c("SYMBOL", "ENTREZID"),
    drop = FALSE
  ]
  converted <- unique(converted)
  converted$list <- list_name
  converted[, c("list", "SYMBOL", "ENTREZID"), drop = FALSE]
}

up_conversion <- symbols_to_entrez(up_symbols, "Up")
down_conversion <- symbols_to_entrez(down_symbols, "Down")
background_conversion <- symbols_to_entrez(
  background_symbols,
  "Background"
)

conversion_report <- rbind(
  up_conversion,
  down_conversion,
  background_conversion
)

up_entrez <- sort(unique(up_conversion$ENTREZID))
down_entrez <- sort(unique(down_conversion$ENTREZID))
background_entrez <- sort(unique(background_conversion$ENTREZID))

# 6. Create mapping and conflict reports --------------------------------------
symbol_mapping_summary <- data.frame(
  list = c("Up", "Down", "Background"),
  input_symbols = c(
    length(up_symbols),
    length(down_symbols),
    length(background_symbols)
  ),
  mapped_symbols = c(
    length(unique(up_conversion$SYMBOL)),
    length(unique(down_conversion$SYMBOL)),
    length(unique(background_conversion$SYMBOL))
  ),
  unique_entrez_ids = c(
    length(up_entrez),
    length(down_entrez),
    length(background_entrez)
  ),
  stringsAsFactors = FALSE
)
symbol_mapping_summary$mapping_rate_percent <- ifelse(
  symbol_mapping_summary$input_symbols == 0,
  NA_real_,
  round(
    100 * symbol_mapping_summary$mapped_symbols /
      symbol_mapping_summary$input_symbols,
    2
  )
)

conflict_report <- deg[
  !is.na(deg$gene_symbol) &
    deg$gene_symbol %in% direction_conflicts,
  c(
    "PROBE_ID", "gene_symbol", "logFC",
    "P.Value", "adj.P.Val", "status"
  ),
  drop = FALSE
]

# 7. Write one-ID-per-line files for DAVID ------------------------------------
writeLines(
  up_symbols,
  file.path(output_dir, "up_gene_symbols.txt")
)
writeLines(
  down_symbols,
  file.path(output_dir, "down_gene_symbols.txt")
)
writeLines(
  all_deg_symbols,
  file.path(output_dir, "all_deg_gene_symbols.txt")
)
writeLines(
  background_symbols,
  file.path(output_dir, "background_gene_symbols.txt")
)
writeLines(
  up_entrez,
  file.path(output_dir, "up_entrez_ids.txt")
)
writeLines(
  down_entrez,
  file.path(output_dir, "down_entrez_ids.txt")
)
writeLines(
  background_entrez,
  file.path(output_dir, "background_entrez_ids.txt")
)

write.csv(
  conversion_report,
  file.path(output_dir, "gene_id_conversion_report.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  symbol_mapping_summary,
  file.path(output_dir, "gene_id_mapping_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  conflict_report,
  file.path(output_dir, "opposite_direction_probe_conflicts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# 8. Record the enrichment-input decision ------------------------------------
input_record <- c(
  "analysis_type: GO/KEGG over-representation analysis preparation",
  "species: Homo sapiens",
  "DEG_source: GSE54568_limma_annotated.csv",
  "direction: Up and Down prepared separately",
  "background: all unique valid gene symbols represented by analyzed probes",
  "multi-probe opposite-direction genes: excluded from Up and Down lists",
  paste0("Up_symbols: ", length(up_symbols)),
  paste0("Down_symbols: ", length(down_symbols)),
  paste0("Conflict_symbols: ", length(direction_conflicts)),
  paste0("Background_symbols: ", length(background_symbols)),
  paste0("Up_ENTREZID: ", length(up_entrez)),
  paste0("Down_ENTREZID: ", length(down_entrez)),
  paste0("Background_ENTREZID: ", length(background_entrez))
)
writeLines(
  input_record,
  file.path(output_dir, "enrichment_input_record.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

print(symbol_mapping_summary)
message("Step 09 complete. Files saved in ", output_dir)
