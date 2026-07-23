# Step 06. Download GEO data and organize sample metadata
# Example: GSE54568 / GPL570

# 1. Install and load packages ------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("GEOquery", "Biobase")
missing_bioc <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

library(GEOquery)
library(Biobase)

# 2. Download GSE -------------------------------------------------------------
options(timeout = 600)
options(download.file.method = "libcurl")

gse_id <- "GSE54568"
target_gpl <- "GPL570"
cache_dir <- "C:/GEO_cache"
output_dir <- "data_processed"

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gset_list <- getGEO(
  gse_id,
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  destdir = cache_dir
)

# 3. Inspect and select GPL ---------------------------------------------------
gpl_summary <- data.frame(
  object = names(gset_list),
  platform = vapply(gset_list, annotation, character(1)),
  probes = vapply(gset_list, nrow, integer(1)),
  samples = vapply(gset_list, ncol, integer(1)),
  stringsAsFactors = FALSE
)
print(gpl_summary)

idx <- which(gpl_summary$platform == target_gpl)
if (length(idx) != 1) {
  stop("분석할 GPL을 하나만 선택하지 못했습니다.")
}

eset <- gset_list[[idx]]

# 4. Separate ExpressionSet components ---------------------------------------
expr_matrix <- exprs(eset)
sample_info <- pData(eset)
feature_info <- fData(eset)

print(dim(expr_matrix))
print(dim(sample_info))
print(dim(feature_info))

# 5. Validate identifiers and order ------------------------------------------
if (!setequal(colnames(expr_matrix), rownames(sample_info))) {
  stop("발현행렬과 sample metadata의 GSM 구성이 다릅니다.")
}

sample_info <- sample_info[
  match(colnames(expr_matrix), rownames(sample_info)),
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(expr_matrix), rownames(sample_info))
)

if (!setequal(rownames(expr_matrix), rownames(feature_info))) {
  stop("발현행렬과 feature annotation의 probe 구성이 다릅니다.")
}

feature_info <- feature_info[
  match(rownames(expr_matrix), rownames(feature_info)),
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(expr_matrix), rownames(feature_info))
)

# 6. Inspect sample metadata --------------------------------------------------
print(colnames(sample_info))

core_cols <- intersect(
  c(
    "geo_accession", "title", "source_name_ch1",
    "organism_ch1", "characteristics_ch1",
    "platform_id"
  ),
  colnames(sample_info)
)
print(head(sample_info[, core_cols, drop = FALSE]))

# 7. Convert characteristics to long format ----------------------------------
char_cols <- grep(
  "^characteristics_ch1",
  colnames(sample_info),
  value = TRUE
)

if (length(char_cols) == 0) {
  warning("characteristics_ch1 열을 찾지 못했습니다.")
  char_long <- data.frame(
    GSM = character(),
    field = character(),
    value = character(),
    raw = character(),
    stringsAsFactors = FALSE
  )
} else {
  char_list <- lapply(seq_len(nrow(sample_info)), function(i) {
    values <- unlist(
      sample_info[i, char_cols, drop = FALSE],
      use.names = FALSE
    )
    values <- trimws(as.character(values))
    values <- values[!is.na(values) & nzchar(values)]

    if (length(values) == 0) {
      return(NULL)
    }

    has_colon <- grepl(":", values, fixed = TRUE)
    fields <- ifelse(
      has_colon,
      trimws(sub(":.*$", "", values)),
      NA_character_
    )
    parsed_values <- ifelse(
      has_colon,
      trimws(sub("^[^:]*:[[:space:]]*", "", values)),
      values
    )

    data.frame(
      GSM = rownames(sample_info)[i],
      field = fields,
      value = parsed_values,
      raw = values,
      stringsAsFactors = FALSE
    )
  })

  char_list <- Filter(Negate(is.null), char_list)
  if (length(char_list) == 0) {
    char_long <- data.frame(
      GSM = character(),
      field = character(),
      value = character(),
      raw = character(),
      stringsAsFactors = FALSE
    )
  } else {
    char_long <- do.call(rbind, char_list)
    rownames(char_long) <- NULL
  }
}

print(head(char_long))
print(sort(unique(stats::na.omit(char_long$field))))

# 8. Create analysis metadata for GSE54568 -----------------------------------
required_meta <- c("title", "source_name_ch1")
missing_meta <- setdiff(required_meta, colnames(sample_info))
if (length(missing_meta) > 0) {
  stop(
    "필수 metadata 열이 없습니다: ",
    paste(missing_meta, collapse = ", ")
  )
}

analysis_meta <- data.frame(
  GSM = rownames(sample_info),
  title = as.character(sample_info$title),
  source = as.character(sample_info$source_name_ch1),
  stringsAsFactors = FALSE
)

analysis_meta$group <- NA_character_
analysis_meta$group[
  grepl("MDD", analysis_meta$title, ignore.case = TRUE)
] <- "Disease"
analysis_meta$group[
  grepl("control", analysis_meta$title, ignore.case = TRUE)
] <- "Control"

print(table(analysis_meta$group, useNA = "ifany"))
print(analysis_meta[, c("GSM", "title", "group")])

if (anyNA(analysis_meta$group)) {
  stop("그룹이 지정되지 않은 sample이 있습니다. metadata를 확인하세요.")
}

analysis_meta$group <- factor(
  analysis_meta$group,
  levels = c("Control", "Disease")
)

stopifnot(
  identical(colnames(expr_matrix), analysis_meta$GSM)
)

# 9. Check repeated subjects when subject ID exists --------------------------
subject_field_names <- c(
  "subject id", "subject_id", "patient id", "patient_id"
)
subject_rows <- char_long[
  !is.na(char_long$field) &
    tolower(char_long$field) %in% subject_field_names,
  c("GSM", "value"),
  drop = FALSE
]

if (nrow(subject_rows) > 0) {
  colnames(subject_rows)[2] <- "subject_id"
  print(table(subject_rows$subject_id))

  repeated_subjects <- subject_rows[
    duplicated(subject_rows$subject_id) |
      duplicated(subject_rows$subject_id, fromLast = TRUE),
    ,
    drop = FALSE
  ]
  print(repeated_subjects)
} else {
  message("subject ID로 판단할 수 있는 characteristics 항목이 없습니다.")
}

# 10. Save original and analysis-ready data ----------------------------------
saveRDS(
  eset,
  file.path(output_dir, "expression_set.rds")
)

expression_export <- data.frame(
  PROBE_ID = rownames(expr_matrix),
  expr_matrix,
  check.names = FALSE
)
write.csv(
  expression_export,
  file.path(output_dir, "expression_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

sample_export <- cbind(
  GSM = rownames(sample_info),
  sample_info
)
write.csv(
  sample_export,
  file.path(output_dir, "sample_metadata_full.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

write.csv(
  analysis_meta,
  file.path(output_dir, "sample_metadata_analysis.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

write.csv(
  char_long,
  file.path(output_dir, "characteristics_long.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

feature_export <- cbind(
  PROBE_ID = rownames(feature_info),
  feature_info
)
write.csv(
  feature_export,
  file.path(output_dir, "feature_annotation.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

write.csv(
  gpl_summary,
  file.path(output_dir, "gpl_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# 11. Read files again and validate -------------------------------------------
eset_check <- readRDS(
  file.path(output_dir, "expression_set.rds")
)
meta_check <- read.csv(
  file.path(output_dir, "sample_metadata_analysis.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  identical(sampleNames(eset_check), meta_check$GSM)
)

print(table(meta_check$group))

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message("Step 06 complete: files saved in ", normalizePath(output_dir))
