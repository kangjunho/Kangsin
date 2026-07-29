# Step 05. GSE69657 Probe ID to Gene Symbol and limma analysis

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

# 2. GSE69657 다운로드 --------------------------------------------------------
gset_list <- getGEO(
  "GSE69657",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  parseCharacteristics = TRUE,
  destdir = "data_raw"
)

print(names(gset_list))
gset <- gset_list[[1]]
print(annotation(gset))
print(gset)

# 3. 발현값, 임상정보, probe annotation --------------------------------------
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

# 4. pData를 직접 확인하고 그룹 열 선택 --------------------------------------
print(colnames(sample_metadata))
View(sample_metadata)

# GSE 설명과 pData를 직접 확인한 뒤 선택한 열입니다.
group_column <- "chemoresponse:ch1"

if (!group_column %in% colnames(sample_metadata)) {
  stop(
    "chemoresponse:ch1 열이 없습니다. ",
    "parseCharacteristics=TRUE로 다시 받고 colnames(sample_metadata)를 확인하세요."
  )
}

print(unique(sample_metadata[[group_column]]))
print(table(sample_metadata[[group_column]], useNA = "ifany"))

chemoresponse <- tolower(
  trimws(as.character(sample_metadata[[group_column]]))
)

group <- factor(
  chemoresponse,
  levels = c("responder", "nonresponder")
)

print(table(group, useNA = "ifany"))

if (anyNA(group)) {
  stop(
    "group에 NA가 있습니다. ",
    "unique(sample_metadata[[group_column]])로 실제 값을 확인하세요."
  )
}

sample_info <- data.frame(
  GSM = rownames(sample_metadata),
  chemoresponse_original = sample_metadata[[group_column]],
  group = group,
  stringsAsFactors = FALSE
)

write.csv(
  sample_info,
  "results/GSE69657_sample_info.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)

# 5. 발현값 범위 확인과 log2 변환 --------------------------------------------
ex <- expression_data

print(
  quantile(
    ex,
    probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
    na.rm = TRUE
  )
)

# GSE69657의 값은 수백~수천 범위이므로 log2 변환합니다.
ex[ex <= 0] <- NA_real_
ex_log2 <- log2(ex)

keep_complete <- rowSums(is.na(ex_log2)) == 0
ex_log2 <- ex_log2[keep_complete, , drop = FALSE]

pdf(
  "results/GSE69657_log2_boxplots.pdf",
  width = 12,
  height = 7
)
par(mfrow = c(2, 1), mar = c(7, 4, 3, 1))
boxplot(
  ex,
  outline = FALSE,
  las = 2,
  main = "Original expression values"
)
boxplot(
  ex_log2,
  outline = FALSE,
  las = 2,
  main = "Log2-transformed expression values"
)
dev.off()

# 6. Probe ID를 Gene symbol로 변환 -------------------------------------------
print(colnames(probe_annotation))

if (!"Gene.symbol" %in% colnames(probe_annotation)) {
  stop(
    "fData에 Gene.symbol 열이 없습니다. ",
    "grep('symbol', colnames(probe_annotation), ignore.case=TRUE, value=TRUE)",
    "로 실제 열 이름을 확인하세요."
  )
}

gene_symbol <- as.character(
  probe_annotation[rownames(ex_log2), "Gene.symbol"]
)

# 여러 symbol이 ///로 연결된 경우 첫 번째 symbol을 대표값으로 사용합니다.
gene_symbol <- trimws(sub("///.*$", "", gene_symbol))

expression_with_symbol <- data.frame(
  Gene.symbol = gene_symbol,
  ex_log2,
  check.names = FALSE
)

# Gene symbol이 없는 probe는 제외합니다.
keep_symbol <-
  !is.na(expression_with_symbol$Gene.symbol) &
  expression_with_symbol$Gene.symbol != "" &
  expression_with_symbol$Gene.symbol != "---"

expression_with_symbol <-
  expression_with_symbol[keep_symbol, , drop = FALSE]

# 같은 Gene symbol에 대응하는 여러 probe의 sample별 평균을 계산합니다.
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

print(dim(ex_log2))
print(dim(gene_expression))
print(head(gene_expression))
View(gene_expression)

stopifnot(
  identical(colnames(gene_expression), rownames(sample_metadata))
)

write.csv(
  gene_expression,
  "results/GSE69657_gene_symbol_expression.csv",
  fileEncoding = "CP949"
)

# 7. Gene symbol 발현행렬로 limma 분석 ---------------------------------------
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

# 8. 분석 환경 기록 ----------------------------------------------------------
capture.output(
  sessionInfo(),
  file = "results/GSE69657_sessionInfo.txt"
)
