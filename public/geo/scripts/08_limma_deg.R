# Step 08. Differential expression analysis with limma
# Example: GSE54568, Disease vs Control

# 1. Install and load limma ---------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
library(limma)

# 2. Load Step 06 and Step 07 outputs ----------------------------------------
data_dir <- "data_processed"
result_dir <- "results/deg"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

expr_probe <- readRDS(
  file.path(data_dir, "expression_probe_processed.rds")
)
sample_info <- read.csv(
  file.path(data_dir, "sample_metadata_analysis.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
probe_map <- read.csv(
  file.path(data_dir, "probe_to_gene_valid.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  is.matrix(expr_probe),
  identical(colnames(expr_probe), sample_info$GSM)
)

# 3. Define and verify groups -------------------------------------------------
sample_info$group <- factor(
  sample_info$group,
  levels = c("Control", "Disease")
)

if (anyNA(sample_info$group)) {
  stop("Control 또는 Disease로 정의되지 않은 sample이 있습니다.")
}

print(table(sample_info$group))
print(sample_info[, c("GSM", "group")])

# 4. Build design and contrast matrices --------------------------------------
design <- model.matrix(
  ~ 0 + group,
  data = sample_info
)
colnames(design) <- sub("^group", "", colnames(design))
rownames(design) <- sample_info$GSM

if (qr(design)$rank < ncol(design)) {
  stop("Design matrix가 full rank가 아닙니다.")
}

contrast_matrix <- makeContrasts(
  Disease_vs_Control = Disease - Control,
  levels = design
)

print(design)
print(colSums(design))
print(contrast_matrix)

# Save the exact model specification.
write.csv(
  cbind(GSM = rownames(design), as.data.frame(design)),
  file.path(result_dir, "GSE54568_design_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  data.frame(
    coefficient = rownames(contrast_matrix),
    contrast_matrix,
    check.names = FALSE
  ),
  file.path(result_dir, "GSE54568_contrast_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# 5. Fit limma model ----------------------------------------------------------
fit_base <- lmFit(expr_probe, design)
fit_contrast <- contrasts.fit(
  fit_base,
  contrast_matrix
)
fit_ebayes <- eBayes(fit_contrast)

# 6. Extract all probe-level results -----------------------------------------
deg_probe <- topTable(
  fit_ebayes,
  coef = "Disease_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg_probe <- data.frame(
  PROBE_ID = rownames(deg_probe),
  deg_probe,
  row.names = NULL,
  check.names = FALSE
)

# 7. Attach gene-symbol annotation without changing result order --------------
map_index <- match(
  deg_probe$PROBE_ID,
  probe_map$probe_id
)
deg_annotated <- deg_probe
deg_annotated$gene_symbol <- probe_map$gene_symbol[map_index]

# 8. Classify significant results --------------------------------------------
fdr_cutoff <- 0.05
lfc_cutoff <- 1

deg_annotated$status <- "Not_significant"
deg_annotated$status[
  !is.na(deg_annotated$adj.P.Val) &
    deg_annotated$adj.P.Val < fdr_cutoff &
    deg_annotated$logFC >= lfc_cutoff
] <- "Up"
deg_annotated$status[
  !is.na(deg_annotated$adj.P.Val) &
    deg_annotated$adj.P.Val < fdr_cutoff &
    deg_annotated$logFC <= -lfc_cutoff
] <- "Down"

deg_sig <- deg_annotated[
  deg_annotated$status != "Not_significant",
  ,
  drop = FALSE
]
deg_up <- deg_annotated[
  deg_annotated$status == "Up",
  ,
  drop = FALSE
]
deg_down <- deg_annotated[
  deg_annotated$status == "Down",
  ,
  drop = FALSE
]

status_summary <- as.data.frame(
  table(deg_annotated$status),
  stringsAsFactors = FALSE
)
colnames(status_summary) <- c("status", "n_probes")
print(status_summary)

# 9. Use treat() to test a minimum effect size -------------------------------
# treat tests whether |logFC| is greater than the specified threshold.
fit_treat <- treat(
  fit_contrast,
  lfc = lfc_cutoff
)
treat_result <- topTreat(
  fit_treat,
  coef = "Disease_vs_Control",
  number = Inf,
  sort.by = "p"
)
treat_result <- data.frame(
  PROBE_ID = rownames(treat_result),
  treat_result,
  row.names = NULL,
  check.names = FALSE
)
treat_result$gene_symbol <- probe_map$gene_symbol[
  match(treat_result$PROBE_ID, probe_map$probe_id)
]

# 10. Sensitivity table for pre-specified thresholds --------------------------
threshold_grid <- expand.grid(
  FDR = c(0.01, 0.05, 0.10),
  abs_logFC = c(0, 0.5, 1),
  stringsAsFactors = FALSE
)

threshold_grid$Up <- mapply(
  function(fdr, lfc) {
    sum(
      deg_annotated$adj.P.Val < fdr &
        deg_annotated$logFC >= lfc,
      na.rm = TRUE
    )
  },
  threshold_grid$FDR,
  threshold_grid$abs_logFC
)
threshold_grid$Down <- mapply(
  function(fdr, lfc) {
    sum(
      deg_annotated$adj.P.Val < fdr &
        deg_annotated$logFC <= -lfc,
      na.rm = TRUE
    )
  },
  threshold_grid$FDR,
  threshold_grid$abs_logFC
)
threshold_grid$Total <- threshold_grid$Up + threshold_grid$Down

# 11. Save results ------------------------------------------------------------
write.csv(
  deg_probe,
  file.path(result_dir, "GSE54568_limma_all_probes.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  deg_annotated,
  file.path(result_dir, "GSE54568_limma_annotated.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  deg_sig,
  file.path(result_dir, "GSE54568_DEG_significant.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  deg_up,
  file.path(result_dir, "GSE54568_upregulated.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  deg_down,
  file.path(result_dir, "GSE54568_downregulated.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  treat_result,
  file.path(result_dir, "GSE54568_treat_logFC1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  threshold_grid,
  file.path(result_dir, "GSE54568_threshold_sensitivity.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
write.csv(
  status_summary,
  file.path(result_dir, "GSE54568_DEG_counts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

saveRDS(
  fit_ebayes,
  file.path(result_dir, "GSE54568_limma_fit.rds")
)

analysis_record <- c(
  "dataset: GSE54568",
  "platform: GPL570",
  "input: expression_probe_processed.rds",
  "comparison: Disease - Control",
  paste0("samples: ", ncol(expr_probe)),
  paste0("probes: ", nrow(expr_probe)),
  "method: limma lmFit + contrasts.fit + eBayes",
  "multiple_testing: Benjamini-Hochberg",
  paste0("FDR_cutoff: ", fdr_cutoff),
  paste0("absolute_logFC_cutoff: ", lfc_cutoff),
  paste0("Up_probes: ", nrow(deg_up)),
  paste0("Down_probes: ", nrow(deg_down))
)
writeLines(
  analysis_record,
  file.path(result_dir, "DEG_analysis_record.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(result_dir, "sessionInfo.txt")
)

message("Step 08 complete. Results saved in ", result_dir)

# -----------------------------------------------------------------------------
# Optional model examples
# Do not run these sections unless the required verified metadata are present.
# -----------------------------------------------------------------------------

if (FALSE) {
  # A. Adjust for age, sex and batch
  required_covariates <- c("age", "sex", "batch")
  stopifnot(all(required_covariates %in% colnames(sample_info)))

  sample_info$age <- as.numeric(sample_info$age)
  sample_info$sex <- factor(sample_info$sex)
  sample_info$batch <- factor(sample_info$batch)

  complete_sample <- complete.cases(
    sample_info[, c("group", required_covariates)]
  )
  meta_cov <- droplevels(sample_info[complete_sample, ])
  expr_cov <- expr_probe[, complete_sample, drop = FALSE]

  design_cov <- model.matrix(
    ~ 0 + group + age + sex + batch,
    data = meta_cov
  )
  colnames(design_cov)[
    grepl("^group", colnames(design_cov))
  ] <- sub(
    "^group",
    "",
    colnames(design_cov)[grepl("^group", colnames(design_cov))]
  )

  if (qr(design_cov)$rank < ncol(design_cov)) {
    stop("공변량 design matrix가 full rank가 아닙니다.")
  }

  contrast_cov <- makeContrasts(
    Disease_vs_Control = Disease - Control,
    levels = design_cov
  )
  fit_cov <- lmFit(expr_cov, design_cov)
  fit_cov <- contrasts.fit(fit_cov, contrast_cov)
  fit_cov <- eBayes(fit_cov)
}

if (FALSE) {
  # B. Repeated measurements from the same subject
  # Use only when each subject ID truly links repeated biological samples.
  stopifnot("subject_id" %in% colnames(sample_info))
  block <- factor(sample_info$subject_id)

  correlation_fit <- duplicateCorrelation(
    expr_probe,
    design,
    block = block
  )
  fit_pair <- lmFit(
    expr_probe,
    design,
    block = block,
    correlation = correlation_fit$consensus
  )
  fit_pair <- contrasts.fit(fit_pair, contrast_matrix)
  fit_pair <- eBayes(fit_pair)
}
