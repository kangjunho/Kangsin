# Step 05. Reproduce a GEO2R analysis in local R
# Example: GSE54568, GPL570, Disease vs Control

# 1. Install missing packages -------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("GEOquery", "limma", "Biobase")
missing_bioc <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

library(GEOquery)
library(limma)
library(Biobase)
library(umap)

# 2. Download GSE -------------------------------------------------------------
options(timeout = 600)
options(download.file.method = "libcurl")

cache_dir <- "C:/GEO_cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

gset_list <- getGEO(
  "GSE54568",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  destdir = cache_dir
)

# 3. Select GPL ---------------------------------------------------------------
print(names(gset_list))

if (length(gset_list) > 1) {
  idx <- grep("GPL570", names(gset_list), fixed = TRUE)
  if (length(idx) != 1) {
    stop("GPL570에 해당하는 ExpressionSet을 하나만 선택하지 못했습니다.")
  }
} else {
  idx <- 1
}

gset <- gset_list[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset), unique = TRUE)

print(annotation(gset))
print(dim(exprs(gset)))
print(dim(pData(gset)))

# 4. Define and verify groups -------------------------------------------------
# One character per sample in the same order as sampleNames(gset).
gsms <- "111111111111111000000000000000"
sml <- strsplit(gsms, split = "")[[1]]

stopifnot(length(sml) == ncol(gset))

groups <- make.names(c("Disease", "Control"))
gs <- factor(sml, levels = c("1", "0"), labels = groups)

print(table(gs, useNA = "ifany"))
group_check <- data.frame(
  GSM = sampleNames(gset),
  group = gs,
  stringsAsFactors = FALSE
)
print(group_check)
write.csv(
  group_check,
  file.path(cache_dir, "GSE54568_group_check.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# 5. Check log2 transformation ------------------------------------------------
ex <- exprs(gset)
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
  exprs(gset) <- log2(ex)
}

keep_complete <- rowSums(is.na(exprs(gset))) == 0
gset <- gset[keep_complete, ]

# 6. Fit limma model ----------------------------------------------------------
gset$group <- gs
design <- model.matrix(~ 0 + group, data = pData(gset))
colnames(design) <- levels(gs)

print(design)
print(colSums(design))

contrast_text <- "Disease-Control"
cont.matrix <- makeContrasts(
  contrasts = contrast_text,
  levels = design
)
print(cont.matrix)

fit <- lmFit(gset, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 7. Export all probe-level results -------------------------------------------
deg_all <- topTable(
  fit2,
  coef = 1,
  adjust.method = "BH",
  sort.by = "P",
  number = Inf
)

print(colnames(deg_all))
print(colnames(deg_all)[duplicated(colnames(deg_all))])

wanted <- c(
  "ID", "logFC", "AveExpr", "t", "P.Value",
  "adj.P.Val", "B", "Gene.symbol",
  "Gene.Symbol", "Gene.title"
)
available <- intersect(wanted, colnames(deg_all))
deg_export <- deg_all[, available, drop = FALSE]

write.csv(
  deg_export,
  file.path(cache_dir, "GSE54568_Disease_vs_Control_all_probes.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# Example DEG threshold: modify before formal analysis.
deg_class <- decideTests(
  fit2,
  adjust.method = "BH",
  p.value = 0.05,
  lfc = 0
)
print(summary(deg_class))

# 8. Save QC plots ------------------------------------------------------------
pdf(
  file.path(cache_dir, "GSE54568_GEO2R_QC.pdf"),
  width = 10,
  height = 7
)

par(mfrow = c(2, 2))
hist(
  deg_all$adj.P.Val,
  col = "grey",
  border = "white",
  xlab = "Adjusted P-value",
  main = "Adjusted P-value distribution"
)
vennDiagram(deg_class, circle.col = "steelblue")
volcanoplot(
  fit2,
  coef = 1,
  main = contrast_text,
  pch = 20,
  highlight = sum(deg_class[, 1] != 0),
  names = rep("+", nrow(fit2))
)
plotMD(
  fit2,
  column = 1,
  status = deg_class[, 1],
  legend = FALSE,
  pch = 20
)
abline(h = 0)

dev.off()

# 9. UMAP ---------------------------------------------------------------------
ex <- exprs(gset)
keep_var <- apply(ex, 1, var, na.rm = TRUE) > 0
ex_umap <- ex[keep_var, , drop = FALSE]

n_nbr <- min(13, ncol(ex_umap) - 1)
if (n_nbr < 2) {
  stop("UMAP을 수행하기에 sample 수가 너무 적습니다.")
}

ump <- umap(
  t(ex_umap),
  n_neighbors = n_nbr,
  random_state = 123
)

pdf(
  file.path(cache_dir, "GSE54568_UMAP.pdf"),
  width = 7,
  height = 6
)
plot(
  ump$layout,
  main = paste0("UMAP, neighbors = ", n_nbr),
  xlab = "UMAP1",
  ylab = "UMAP2",
  col = as.integer(gs),
  pch = 20,
  cex = 1.5
)
legend(
  "topright",
  legend = levels(gs),
  col = seq_along(levels(gs)),
  pch = 20,
  title = "Group"
)
dev.off()

# 10. Reproducibility ---------------------------------------------------------
capture.output(
  sessionInfo(),
  file = file.path(cache_dir, "GSE54568_sessionInfo.txt")
)
