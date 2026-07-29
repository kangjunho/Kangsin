# Step 05. GSE69657 직접 전처리와 gene-level limma 분석

# 1. 패키지와 폴더 -----------------------------------------------------------
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

# 2. GSE69657 Series Matrix ---------------------------------------------------
gset_list <- getGEO(
  "GSE69657",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  destdir = "data_raw"
)

print(names(gset_list))
gset <- gset_list[[1]]
print(annotation(gset))
print(gset)

# 3. exprs, pData, fData ------------------------------------------------------
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
  fileEncoding = "CP949"
)
write.csv(
  sample_metadata,
  "data_raw/GSE69657_pData.csv",
  fileEncoding = "CP949"
)
write.csv(
  probe_annotation,
  "data_raw/GSE69657_fData.csv",
  fileEncoding = "CP949"
)

# 4. chemoresponse 그룹 -------------------------------------------------------
# GEO 화면의 chemoresponse:ch1 정보는 GEOquery에서
# characteristics_ch1.3 같은 이름으로 표시될 수 있습니다.
response_candidates <- colnames(sample_metadata)[
  vapply(
    sample_metadata,
    function(x) {
      any(grepl(
        "^chemoresponse:",
        as.character(x),
        ignore.case = TRUE
      ))
    },
    logical(1)
  )
]

print(response_candidates)

if (length(response_candidates) != 1) {
  stop(
    "chemoresponse 후보 열이 하나가 아닙니다. ",
    "colnames(sample_metadata)와 unique() 값을 확인하세요."
  )
}

response_col <- response_candidates[1]
print(unique(sample_metadata[[response_col]]))

chemoresponse <- trimws(
  sub(
    "^chemoresponse:\\s*",
    "",
    sample_metadata[[response_col]],
    ignore.case = TRUE
  )
)

group <- factor(
  chemoresponse,
  levels = c("responder", "nonresponder")
)

print(table(group, useNA = "ifany"))

if (anyNA(group)) {
  stop("group에 NA가 있습니다. chemoresponse의 실제 철자를 확인하세요.")
}

sample_info <- data.frame(
  GSM = rownames(sample_metadata),
  chemoresponse_original = sample_metadata[[response_col]],
  group = group,
  stringsAsFactors = FALSE
)

write.csv(
  sample_info,
  "results/GSE69657_sample_info.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)

# 5. 30개 Series Matrix 전처리 -----------------------------------------------
ex <- expression_data

print(
  quantile(
    ex,
    probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
    na.rm = TRUE
  )
)

ex[ex <= 0] <- NA_real_
ex_log2 <- log2(ex)

keep_complete <- rowSums(is.na(ex_log2)) == 0
ex_log2 <- ex_log2[keep_complete, , drop = FALSE]

ex_normalized <- normalizeBetweenArrays(
  ex_log2,
  method = "quantile"
)

pdf(
  "results/GSE69657_normalization_boxplots.pdf",
  width = 12,
  height = 7
)
par(mfrow = c(2, 1), mar = c(7, 4, 3, 1))
boxplot(
  ex_log2,
  outline = FALSE,
  las = 2,
  main = "Before quantile normalization"
)
boxplot(
  ex_normalized,
  outline = FALSE,
  las = 2,
  main = "After quantile normalization"
)
dev.off()

# 6. Probe ID를 gene symbol로 평균 통합 --------------------------------------
if (!"Gene.symbol" %in% colnames(probe_annotation)) {
  stop(
    "fData에 Gene.symbol 열이 없습니다. ",
    "grep('symbol', colnames(probe_annotation), ignore.case=TRUE, value=TRUE)",
    "로 실제 열 이름을 확인하세요."
  )
}

gene_symbol <- as.character(
  probe_annotation[rownames(ex_normalized), "Gene.symbol"]
)
gene_symbol <- trimws(sub("///.*$", "", gene_symbol))

expression_with_symbol <- data.frame(
  Gene.symbol = gene_symbol,
  ex_normalized,
  check.names = FALSE
)

keep_symbol <-
  !is.na(expression_with_symbol$Gene.symbol) &
  expression_with_symbol$Gene.symbol != "" &
  expression_with_symbol$Gene.symbol != "---"

expression_with_symbol <-
  expression_with_symbol[keep_symbol, , drop = FALSE]

gene_expression_df <- aggregate(
  . ~ Gene.symbol,
  data = expression_with_symbol,
  FUN = mean,
  na.rm = TRUE
)

rownames(gene_expression_df) <- gene_expression_df$Gene.symbol
gene_expression <- as.matrix(
  gene_expression_df[, -1, drop = FALSE]
)
storage.mode(gene_expression) <- "numeric"

print(dim(ex_normalized))
print(dim(gene_expression))

write.csv(
  gene_expression,
  "results/GSE69657_gene_expression_normalized.csv",
  fileEncoding = "CP949"
)

# 7. Gene-level limma ---------------------------------------------------------
stopifnot(
  identical(colnames(gene_expression), rownames(sample_metadata)),
  length(group) == ncol(gene_expression)
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Responder", "Nonresponder")
rownames(design) <- colnames(gene_expression)

print(design)
print(colSums(design))

contrast_matrix <- makeContrasts(
  Responder_vs_Nonresponder =
    Responder - Nonresponder,
  levels = design
)
print(contrast_matrix)

fit <- lmFit(gene_expression, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(
  fit2,
  coef = "Responder_vs_Nonresponder",
  adjust.method = "BH",
  sort.by = "P",
  number = Inf
)
deg_results$Gene.symbol <- rownames(deg_results)

up_degs <- deg_results[
  deg_results$P.Value < 0.05 &
    deg_results$logFC > 1,
  ,
  drop = FALSE
]

down_degs <- deg_results[
  deg_results$P.Value < 0.05 &
    deg_results$logFC < -1,
  ,
  drop = FALSE
]

write.csv(
  deg_results,
  "results/GSE69657_all_genes.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)
write.csv(
  up_degs,
  "results/GSE69657_up_degs.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)
write.csv(
  down_degs,
  "results/GSE69657_down_degs.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)

# 8. 선택 실습: 16개 raw CEL에 RMA -------------------------------------------
# 주의: GSE69657의 30개 sample 중 CEL 파일은 16개만 제공됩니다.
# 아래 분석은 30개 Series Matrix 분석과 별도의 부분집합 분석입니다.
#
# BiocManager::install(c("affy", "hgu133plus2cdf"))
# install.packages("R.utils")
# library(affy)
# library(R.utils)
#
# getGEOSuppFiles(
#   "GSE69657",
#   makeDirectory = TRUE,
#   baseDir = "data_raw"
# )
#
# untar(
#   "data_raw/GSE69657/GSE69657_RAW.tar",
#   exdir = "data_raw/GSE69657/CEL"
# )
#
# cel_gz <- list.files(
#   "data_raw/GSE69657/CEL",
#   pattern = "\\.CEL\\.gz$",
#   full.names = TRUE,
#   ignore.case = TRUE
# )
#
# vapply(
#   cel_gz,
#   gunzip,
#   FUN.VALUE = character(1),
#   remove = FALSE,
#   overwrite = TRUE
# )
#
# cel_files <- list.files(
#   "data_raw/GSE69657/CEL",
#   pattern = "\\.CEL$",
#   full.names = TRUE,
#   ignore.case = TRUE
# )
# stopifnot(length(cel_files) == 16)
#
# raw_affy <- ReadAffy(filenames = cel_files)
# rma_eset <- rma(raw_affy)
# rma_expression <- exprs(rma_eset)
# print(dim(rma_expression))
#
# rma_gene_symbol <- as.character(
#   probe_annotation[rownames(rma_expression), "Gene.symbol"]
# )
# rma_gene_symbol <- trimws(sub("///.*$", "", rma_gene_symbol))
#
# rma_with_symbol <- data.frame(
#   Gene.symbol = rma_gene_symbol,
#   rma_expression,
#   check.names = FALSE
# )
# rma_with_symbol <- rma_with_symbol[
#   !is.na(rma_with_symbol$Gene.symbol) &
#     rma_with_symbol$Gene.symbol != "" &
#     rma_with_symbol$Gene.symbol != "---",
#   ,
#   drop = FALSE
# ]
#
# rma_gene_df <- aggregate(
#   . ~ Gene.symbol,
#   data = rma_with_symbol,
#   FUN = mean,
#   na.rm = TRUE
# )
# rownames(rma_gene_df) <- rma_gene_df$Gene.symbol
# rma_gene_expression <- as.matrix(
#   rma_gene_df[, -1, drop = FALSE]
# )

# 9. 분석 환경 기록 ----------------------------------------------------------
capture.output(
  sessionInfo(),
  file = "results/GSE69657_sessionInfo.txt"
)
