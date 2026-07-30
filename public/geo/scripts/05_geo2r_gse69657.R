# Step 05. GSE69657 직접 분석과 GEO2R 결과 비교

# Step 01에서 설치한 패키지를 불러옵니다.
library(GEOquery)
library(Biobase)
library(limma)

# 1. GSE69657 다운로드 --------------------------------------------------------
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

# 2. 발현값, 임상정보, probe annotation --------------------------------------
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

# 3. pData를 직접 확인하고 그룹 열 선택 --------------------------------------
print(colnames(sample_metadata))
View(sample_metadata)

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

# 이 데이터의 실제 표기는 "noresponder"입니다.
group <- factor(
  chemoresponse,
  levels = c("responder", "noresponder")
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

# 4. 발현값 범위 확인과 log2 변환 --------------------------------------------
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

pdf("results/GSE69657_log2_boxplots.pdf", width = 12, height = 7)
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

# 5. GEO2R 방식의 probe-level limma 결과 -------------------------------------
probe_design <- model.matrix(~ 0 + group)
colnames(probe_design) <- c("Responder", "Noresponder")
rownames(probe_design) <- colnames(ex_log2)

probe_contrast <- makeContrasts(
  Responder_vs_Noresponder = Responder - Noresponder,
  levels = probe_design
)

probe_fit <- lmFit(ex_log2, probe_design)
probe_fit2 <- contrasts.fit(probe_fit, probe_contrast)
probe_fit2 <- eBayes(probe_fit2, proportion = 0.01)

geo2r_results <- topTable(
  probe_fit2,
  coef = "Responder_vs_Noresponder",
  adjust.method = "BH",
  sort.by = "B",
  number = Inf
)
geo2r_results$ID <- rownames(geo2r_results)

# 실행 결과에서 확인된 실제 annotation 열 이름은 "Gene symbol"입니다.
symbol_column <- "Gene symbol"
if (!symbol_column %in% colnames(probe_annotation)) {
  stop(
    "fData에 'Gene symbol' 열이 없습니다. ",
    "grep('symbol', colnames(probe_annotation), ignore.case=TRUE, value=TRUE)",
    "로 실제 열 이름을 확인하세요."
  )
}

geo2r_results$Gene.symbol <- as.character(
  probe_annotation[geo2r_results$ID, symbol_column]
)
geo2r_results$Gene.symbol <- trimws(
  sub("///.*$", "", geo2r_results$Gene.symbol)
)

needed_cols <- c(
  "ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "Gene.symbol"
)
geo2r_results <- geo2r_results[, needed_cols, drop = FALSE]

write.csv(
  geo2r_results,
  "results/GSE69657_GEO2R_probe_results.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)

# 6. Probe ID를 Gene symbol로 변환 -------------------------------------------
gene_symbol <- as.character(
  probe_annotation[rownames(ex_log2), symbol_column]
)

# 여러 symbol이 ///로 연결된 경우 첫 번째 symbol을 대표값으로 사용합니다.
gene_symbol <- trimws(sub("///.*$", "", gene_symbol))

expression_with_symbol <- data.frame(
  Gene.symbol = gene_symbol,
  ex_log2,
  check.names = FALSE
)

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
gene_design <- model.matrix(~ 0 + group)
colnames(gene_design) <- c("Responder", "Noresponder")
rownames(gene_design) <- colnames(gene_expression)

print(gene_design)
print(colSums(gene_design))

gene_contrast <- makeContrasts(
  Responder_vs_Noresponder = Responder - Noresponder,
  levels = gene_design
)

gene_fit <- lmFit(gene_expression, gene_design)
gene_fit2 <- contrasts.fit(gene_fit, gene_contrast)
gene_fit2 <- eBayes(gene_fit2)

deg_results <- topTable(
  gene_fit2,
  coef = "Responder_vs_Noresponder",
  adjust.method = "BH",
  sort.by = "P",
  number = Inf
)
deg_results$Gene.symbol <- rownames(deg_results)

# 실습용 기준: nominal P-value < 0.05, |logFC| > 1
up_degs <- deg_results[
  !is.na(deg_results$P.Value) &
    !is.na(deg_results$logFC) &
    deg_results$P.Value < 0.05 &
    deg_results$logFC > 1,
  ,
  drop = FALSE
]

down_degs <- deg_results[
  !is.na(deg_results$P.Value) &
    !is.na(deg_results$logFC) &
    deg_results$P.Value < 0.05 &
    deg_results$logFC < -1,
  ,
  drop = FALSE
]

write.csv(
  deg_results,
  "results/GSE69657_gene_level_all_results.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)
write.csv(
  up_degs,
  "results/GSE69657_gene_level_up_degs.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)
write.csv(
  down_degs,
  "results/GSE69657_gene_level_down_degs.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)

# 8. GEO2R와 직접 분석 DEG 결과 비교 ----------------------------------------
# GEO2R 결과는 probe 단위이고 직접 분석은 gene 단위입니다.
# 같은 유전자에 여러 probe가 있으면 P.Value가 가장 작은 probe를
# 그 유전자의 대표 probe로 선택한 뒤 비교합니다.
geo2r_gene <- geo2r_results[
  !is.na(geo2r_results$Gene.symbol) &
    geo2r_results$Gene.symbol != "" &
    geo2r_results$Gene.symbol != "---",
  ,
  drop = FALSE
]
geo2r_gene <- geo2r_gene[
  order(geo2r_gene$Gene.symbol, geo2r_gene$P.Value),
  ,
  drop = FALSE
]
geo2r_gene <- geo2r_gene[
  !duplicated(geo2r_gene$Gene.symbol),
  ,
  drop = FALSE
]

geo2r_gene$GEO2R.class <- "Not_DEG"
geo2r_gene$GEO2R.class[
  geo2r_gene$P.Value < 0.05 & geo2r_gene$logFC > 1
] <- "Up"
geo2r_gene$GEO2R.class[
  geo2r_gene$P.Value < 0.05 & geo2r_gene$logFC < -1
] <- "Down"

direct_gene <- deg_results
direct_gene$Direct.class <- "Not_DEG"
direct_gene$Direct.class[
  direct_gene$P.Value < 0.05 & direct_gene$logFC > 1
] <- "Up"
direct_gene$Direct.class[
  direct_gene$P.Value < 0.05 & direct_gene$logFC < -1
] <- "Down"

comparison <- merge(
  geo2r_gene[, c(
    "Gene.symbol", "ID", "logFC", "P.Value", "adj.P.Val", "GEO2R.class"
  )],
  direct_gene[, c(
    "Gene.symbol", "logFC", "P.Value", "adj.P.Val", "Direct.class"
  )],
  by = "Gene.symbol",
  suffixes = c(".GEO2R", ".Direct")
)

comparison$Direction.same <- with(
  comparison,
  sign(logFC.GEO2R) == sign(logFC.Direct)
)

print(table(comparison$GEO2R.class, comparison$Direct.class))
print(cor(
  comparison$logFC.GEO2R,
  comparison$logFC.Direct,
  use = "complete.obs",
  method = "pearson"
))
print(mean(comparison$Direction.same, na.rm = TRUE))

geo2r_up <- comparison$Gene.symbol[comparison$GEO2R.class == "Up"]
geo2r_down <- comparison$Gene.symbol[comparison$GEO2R.class == "Down"]
direct_up <- comparison$Gene.symbol[comparison$Direct.class == "Up"]
direct_down <- comparison$Gene.symbol[comparison$Direct.class == "Down"]

jaccard <- function(x, y) {
  union_set <- union(x, y)
  if (length(union_set) == 0) return(NA_real_)
  length(intersect(x, y)) / length(union_set)
}

comparison_summary <- data.frame(
  Metric = c(
    "Compared genes",
    "logFC Pearson correlation",
    "Same-direction proportion",
    "GEO2R Up genes",
    "Direct Up genes",
    "Common Up genes",
    "Up Jaccard index",
    "GEO2R Down genes",
    "Direct Down genes",
    "Common Down genes",
    "Down Jaccard index"
  ),
  Value = c(
    nrow(comparison),
    cor(
      comparison$logFC.GEO2R,
      comparison$logFC.Direct,
      use = "complete.obs",
      method = "pearson"
    ),
    mean(comparison$Direction.same, na.rm = TRUE),
    length(geo2r_up),
    length(direct_up),
    length(intersect(geo2r_up, direct_up)),
    jaccard(geo2r_up, direct_up),
    length(geo2r_down),
    length(direct_down),
    length(intersect(geo2r_down, direct_down)),
    jaccard(geo2r_down, direct_down)
  )
)
print(comparison_summary)

write.csv(
  comparison,
  "results/GSE69657_GEO2R_vs_gene_level.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)
write.csv(
  comparison_summary,
  "results/GSE69657_GEO2R_vs_gene_level_summary.csv",
  row.names = FALSE,
  fileEncoding = "CP949"
)

pdf("results/GSE69657_GEO2R_vs_gene_level_logFC.pdf")
plot(
  comparison$logFC.GEO2R,
  comparison$logFC.Direct,
  pch = 16,
  cex = 0.5,
  col = rgb(0.1, 0.4, 0.8, 0.35),
  xlab = "GEO2R representative-probe logFC",
  ylab = "Direct gene-level logFC",
  main = "GSE69657: GEO2R vs direct analysis"
)
abline(h = 0, v = 0, col = "grey70", lty = 2)
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()
