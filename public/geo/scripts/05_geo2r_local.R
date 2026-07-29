# Step 05. Direct analysis of GSE69657 and comparison with GEO2R

# 1. Install and load packages ------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

needed <- c("GEOquery", "Biobase", "limma")
missing <- needed[
  !vapply(needed, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0) {
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

library(GEOquery)
library(Biobase)
library(limma)

dir.create("data_raw", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
options(timeout = 600)

# 2. Download and select the ExpressionSet -----------------------------------
gset_list <- getGEO(
  "GSE69657",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  destdir = "data_raw"
)

print(names(gset_list))

if (length(gset_list) > 1) {
  idx <- grep("GPL570", names(gset_list), fixed = TRUE)
  if (length(idx) != 1) {
    stop("GPL570 ExpressionSet을 하나만 선택하지 못했습니다.")
  }
} else {
  idx <- 1
}

gset <- gset_list[[idx]]
print(annotation(gset))
print(gset)

# 3. Extract expression, sample metadata and probe annotation ----------------
expression_data <- exprs(gset)
sample_metadata <- pData(gset)
probe_annotation <- fData(gset)

print(dim(expression_data))
print(dim(sample_metadata))
print(dim(probe_annotation))

stopifnot(
  identical(colnames(expression_data), rownames(sample_metadata)),
  identical(rownames(expression_data), rownames(probe_annotation))
)

write.csv(
  expression_data,
  "data_raw/GSE69657_expression.csv",
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  sample_metadata,
  "data_raw/GSE69657_pData.csv",
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  probe_annotation,
  "data_raw/GSE69657_fData.csv",
  fileEncoding = "UTF-8-BOM"
)

# 4. Find the response column and define groups ------------------------------
print(colnames(sample_metadata))

candidate_columns <- grep(
  "response|characteristics|title|source",
  colnames(sample_metadata),
  ignore.case = TRUE,
  value = TRUE
)
print(candidate_columns)

response_candidates <- colnames(sample_metadata)[
  vapply(
    sample_metadata,
    function(x) {
      any(grepl(
        "responder",
        as.character(x),
        ignore.case = TRUE
      ))
    },
    logical(1)
  )
]

if (length(response_candidates) == 0) {
  stop(
    paste(
      "값에 responder가 포함된 열을 자동으로 찾지 못했습니다.",
      "candidate_columns의 characteristics 열을 직접 확인한 뒤",
      "response_col에 정확한 열 이름을 입력하세요."
    )
  )
}

if (length(response_candidates) > 1) {
  print(response_candidates)
  stop(
    paste(
      "response 후보 열이 두 개 이상입니다.",
      "unique(sample_metadata[[열 이름]])로 값을 확인한 뒤",
      "response_col을 하나만 직접 선택하세요."
    )
  )
}

response_col <- response_candidates[1]
print(response_col)
print(unique(sample_metadata[[response_col]]))

response_text <- tolower(
  trimws(as.character(sample_metadata[[response_col]]))
)

# Check non-responder first because it contains the substring "responder".
group_text <- ifelse(
  grepl("non[-_ ]?responder", response_text),
  "Non_responder",
  ifelse(grepl("responder", response_text), "Responder", NA)
)

group <- factor(
  group_text,
  levels = c("Responder", "Non_responder")
)

print(table(group, useNA = "ifany"))

if (anyNA(group)) {
  stop(
    paste(
      "그룹으로 변환하지 못한 값이 있습니다.",
      "unique(sample_metadata[[response_col]])를 확인하고",
      "group_text 변환 규칙을 실제 표기법에 맞게 수정하세요."
    )
  )
}

group_check <- data.frame(
  GSM = rownames(sample_metadata),
  original_value = sample_metadata[[response_col]],
  group = group,
  stringsAsFactors = FALSE
)
print(group_check)

write.csv(
  group_check,
  "results/GSE69657_group_check.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# 5. Reproduce GEO2R log2 decision -------------------------------------------
ex <- expression_data
qx <- quantile(
  ex,
  probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
  na.rm = TRUE
)
print(qx)

log_needed <- (qx[5] > 100) ||
  ((qx[6] - qx[1] > 50) && qx[2] > 0)

if (log_needed) {
  ex[ex <= 0] <- NA_real_
  ex <- log2(ex)
}

keep <- rowSums(is.na(ex)) == 0
ex <- ex[keep, , drop = FALSE]
probe_annotation <- probe_annotation[rownames(ex), , drop = FALSE]

# 6. Fit the limma model ------------------------------------------------------
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(ex)

print(design)
print(colSums(design))

contrast_matrix <- makeContrasts(
  Responder_vs_Non_responder =
    Responder - Non_responder,
  levels = design
)
print(contrast_matrix)

fit <- lmFit(ex, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, proportion = 0.01)

local_results <- topTable(
  fit2,
  coef = "Responder_vs_Non_responder",
  adjust.method = "BH",
  sort.by = "B",
  number = Inf
)
local_results$ID <- rownames(local_results)

# 7. Add fData annotation -----------------------------------------------------
annotation_match <- probe_annotation[
  match(local_results$ID, rownames(probe_annotation)),
  ,
  drop = FALSE
]

stopifnot(
  identical(local_results$ID, rownames(annotation_match))
)

local_annotated <- cbind(
  local_results,
  annotation_match
)

print(
  grep(
    "gene.*symbol|symbol",
    colnames(local_annotated),
    ignore.case = TRUE,
    value = TRUE
  )
)

write.csv(
  local_annotated,
  "results/GSE69657_local_limma_all_probes.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# 8. Compare with the GEO2R full table ---------------------------------------
# Download "Download full table" from GEO2R and save it as:
geo2r_file <- "data_raw/GSE69657_GEO2R_results.tsv"

if (!file.exists(geo2r_file)) {
  message(
    "Local limma analysis is complete. ",
    "Download the GEO2R full table to ",
    geo2r_file,
    " and rerun the comparison section."
  )
} else {
  geo2r_results <- read.delim(
    geo2r_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_columns <- c("ID", "logFC", "P.Value", "adj.P.Val")
  missing_local <- setdiff(required_columns, colnames(local_results))
  missing_geo2r <- setdiff(required_columns, colnames(geo2r_results))

  if (length(missing_local) > 0) {
    stop(
      "Local result is missing: ",
      paste(missing_local, collapse = ", ")
    )
  }
  if (length(missing_geo2r) > 0) {
    stop(
      "GEO2R result is missing: ",
      paste(missing_geo2r, collapse = ", "),
      ". Check colnames(geo2r_results)."
    )
  }

  local_compare <- local_results[, required_columns]
  colnames(local_compare) <- c(
    "ID", "logFC_local", "P.Value_local", "adj.P.Val_local"
  )

  geo2r_compare <- geo2r_results[, required_columns]
  colnames(geo2r_compare) <- c(
    "ID", "logFC_GEO2R", "P.Value_GEO2R", "adj.P.Val_GEO2R"
  )

  comparison <- merge(
    local_compare,
    geo2r_compare,
    by = "ID",
    all = FALSE
  )

  comparison$logFC_difference <-
    comparison$logFC_local - comparison$logFC_GEO2R

  comparison$same_direction <-
    sign(comparison$logFC_local) ==
    sign(comparison$logFC_GEO2R)

  comparison_summary <- data.frame(
    matched_probes = nrow(comparison),
    logFC_correlation = cor(
      comparison$logFC_local,
      comparison$logFC_GEO2R,
      use = "complete.obs"
    ),
    maximum_absolute_logFC_difference = max(
      abs(comparison$logFC_difference),
      na.rm = TRUE
    ),
    same_direction_proportion = mean(
      comparison$same_direction,
      na.rm = TRUE
    )
  )
  print(comparison_summary)

  write.csv(
    comparison,
    "results/GSE69657_local_vs_GEO2R.csv",
    row.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  write.csv(
    comparison_summary,
    "results/GSE69657_comparison_summary.csv",
    row.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )

  png(
    "results/GSE69657_local_vs_GEO2R_logFC.png",
    width = 1400,
    height = 1200,
    res = 180
  )
  plot(
    comparison$logFC_GEO2R,
    comparison$logFC_local,
    pch = 20,
    col = rgb(0.1, 0.4, 0.8, 0.25),
    xlab = "GEO2R logFC",
    ylab = "Local limma logFC",
    main = "Local limma vs GEO2R"
  )
  abline(a = 0, b = 1, col = "red", lwd = 2)
  dev.off()
}

# 9. Reproducibility ----------------------------------------------------------
capture.output(
  sessionInfo(),
  file = "results/GSE69657_sessionInfo.txt"
)
